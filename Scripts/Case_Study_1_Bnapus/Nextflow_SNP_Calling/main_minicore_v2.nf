
nextflow.enable.dsl=2

/*
 * =========================
 * BASE DIRECTORY
 * =========================
 */

def baseDir = workflow.projectDir

/*
 * =========================
 * INPUTS
 * =========================
 */

params.reads       = "${baseDir}/data/*_{f1,r2}.fq.gz"
params.reference   = "${baseDir}/reference/DarmorV10_Chromosomes_Only.fa"
params.index	   = "${baseDir}/reference/DarmorV10_Chromosomes_Only.fa"
params.known_sites = "${baseDir}/reference/known_sites_DarmorV10.sorted.vcf.gz"
params.intervals   = "${baseDir}/reference/intervals/*.interval_list"

/*
 * =========================
 * OUTPUT BASE
 * =========================
 */

params.outdir      = "${baseDir}/results"

/*
 * =========================
 * BASE OUTPUT STRUCTURE
 * =========================
 */

params.qc           = "${params.outdir}/qc/fastp"
params.alignments   = "${params.outdir}/alignments"
params.sorted	    = "${params.outdir}/alignments/sorted"
params.dedup	    = "${params.outdir}/alignments/marked_duplicates"
params.recal	    = "${params.outdir}/alignments/recalibrated"
params.gvcf         = "${params.outdir}/gvcf/per_sample"
params.genomics_db  = "${params.outdir}/genomics_db/shards"
params.joint        = "${params.outdir}/joint_genotyping"
params.merged       = "${params.outdir}/joint_genotyping/merged"
params.logs         = "${params.outdir}/logs"

/*
 * =========================
 * GENERAL SETTINGS
 * =========================
 */
 
params.vcf_prefix   = "Samples_Bna_Minicore_1-135"
params.db_mode = "create" // or "update"


workflow {

    /*
     * =========================
     * INPUT CHANNELS
     * =========================
     */
	 
	reads_ch = Channel.fromFilePairs(params.reads)
		.map { id, reads ->
			tuple(id, reads[0], reads[1])
		}

	intervals_ch = Channel.fromPath(params.intervals)
		.map { file -> 
			def id = file.simpleName.tokenize('-')[0]
			tuple(id, file)
		}
		
    /*
     * =========================
     * QC (FASTP)
     * =========================
     */
	 
    clean_reads = fastp(reads_ch)
        .map { id, r1, r2, ->
            tuple(id, r1, r2)
        }
	
    /*
     * =========================
     * ALIGNMENT
     * =========================
     */

	ref_index = file(params.index + ".*")
	mapped = mapping(clean_reads, ref_index)
	
    /*
     * =========================
     * MARK DUPLICATES
     * =========================
     */

	marked_bam = mark_duplicates(mapped)
		.map { id, bam, bai, metrics ->
			tuple(id, bam, bai)
		}

    /*
     * =========================
     * PREPARE REFERENCE FILES
     * =========================
     */
	
	known_sites_index = files(params.known_sites + "*")

     /*
      * =========================
      * BQSR
      * =========================
      */

	recal_bam = bqsr(marked_bam, ref_index, known_sites_index)
	
	/*
     * =========================
     * CLEAN BQSR OUTPUT
     * =========================
     */
	
	recal_bam_clean = recal_bam
		.map { id, bam, bai, before, after, csv ->
			tuple(id, bam, bai)
		}
	
	/*
     * =========================
     * SCATTER ACROSS INTERVALS
     * =========================
     */
	
	recal_with_intervals = recal_bam_clean.combine(intervals_ch)
	
	/*
     * =========================
     * HAPLOTYPE CALLER
     * =========================
     */
	
	gvcf = haplotype_caller(recal_with_intervals, ref_index)
	
	/*
	 * =========================
	 * GROUP GVCFs BY INTERVAL
	 * =========================
	 */
	
	gvcf_by_interval = gvcf	
		.map { id, interval_id, gvcf_file, tbi ->
			tuple(interval_id, gvcf_file, tbi)
		}
		.groupTuple()
		.join(intervals_ch)
		.map { interval_id, gvcfs, tbis, interval_file ->
			tuple(interval_id, gvcfs, tbis, interval_file)
		}
	
	/*
     * =========================
     * GENOMICSDB IMPORT
     * =========================
     */
	
	genomicsdb_ch = genomicsdb_import(gvcf_by_interval)
	
	/*
     * =========================
     * REATTACH INTERVAL METADATA
     * =========================
     */
	
	genomicsdb_ready =
		genomicsdb_ch
			.map { interval_id, db ->
				tuple(interval_id, db)
			}
			.join(intervals_ch)
			.map { interval_id, db, interval_file ->
				tuple(interval_id, db, interval_file)
			}
		
	
	/*
     * =========================
     * JOINT CALLING
     * =========================
     */
	
	genotyped = genotype_gvcfs(genomicsdb_ready)
	
	/*
	 * =========================
	 * FINAL MERGE (ORDERED IN SHELL)
	 * =========================
	 */

	vcf_files_ch = genotyped
		.map { interval_id, vcf, tbi -> tuple(vcf, tbi) }
		.collect()
	
	final_vcf = concat_vcfs(vcf_files_ch)
	
}

process fastp {

	tag "$id"
	
	container 'oras://community.wave.seqera.io/library/bcftools_bwa_fastp_gatk4_samtools:0fe683a5007bf05c' 
		
	input:
	tuple val(id), path(r1), path(r2)
	
	output:
	tuple val(id),
		path("${id}_f1.clean.fq.gz"),
		path("${id}_r2.clean.fq.gz")
	
	script:
	"""
	fastp \
		-i ${r1} \
		-I ${r2} \
		-o ${id}_f1.clean.fq.gz \
		-O ${id}_r2.clean.fq.gz \
		--detect_adapter_for_pe \
		--correction \
		--thread 8

	echo ""
	echo "Processed: ${id}"
	echo ""
	echo "In1: ${r1}"
	echo "In2: ${r2}"
	echo ""
	echo "Out1: ${id}_f1.clean.fq.gz"
	echo "Out2: ${id}_r2.clean.fq.gz"
	echo ""
	"""
}

process mapping {
	
	tag "$id"
	
	container 'oras://community.wave.seqera.io/library/bcftools_bwa_fastp_gatk4_samtools:0fe683a5007bf05c'  
		
	input:
	tuple val(id), path(r1), path(r2)
	path(ref_index_files)
	
	output:
	tuple val(id),
		path("${id}.sorted.bam"),
		path("${id}.sorted.bam.bai")
	
	script:
	"""
	
	echo "Mapping sample ${id}"
	
    bwa mem \
	-t 32 \
        -R "@RG\\tID:${id}\\tSM:${id}\\tPL:ILLUMINA\\tLB:${id}\\tPU:${id}" \
        ${params.index} \
        ${r1} \
        ${r2} | \
    samtools sort \
        -@ 16 \
        -o ${id}.sorted.bam

    samtools index ${id}.sorted.bam
	
	echo "Finished mapping: ${id}"
	"""
}
	
process mark_duplicates {
	
	tag "$id"
	
	container 'oras://community.wave.seqera.io/library/bcftools_bwa_fastp_gatk4_samtools:0fe683a5007bf05c'  
	
	publishDir "${params.dedup}", mode: 'copy'
	
	input:
	tuple val(id), path(bam), path(bai)
	
	output:
	tuple val(id),
		path("${id}.marked.bam"),
		path("${id}.marked.bam.bai"),
		path("${id}.dup_metrics.txt")
	
	script:
	"""
	echo "Marking duplicates for ${id}"
	
	mkdir -p tmp	

	gatk --java-options "-Xmx48g" MarkDuplicates \
		-I ${bam} \
		-O ${id}.marked.bam \
		-M ${id}.dup_metrics.txt \
		--CREATE_INDEX false \
		--VALIDATION_STRINGENCY SILENT \
		--TMP_DIR tmp

	gatk BuildBamIndex \
		-I ${id}.marked.bam \
		-O ${id}.marked.bam.bai
	
	echo "Finished duplicate marking for ${id}"
	"""
}

process bqsr {
	
	tag "$id"
	
	container 'oras://community.wave.seqera.io/library/bcftools_bwa_fastp_gatk4_samtools:0fe683a5007bf05c'  
	
	publishDir "${params.recal}", mode: 'copy'
	
	input:
	tuple val(id), path(marked_bam), path(bai)
	path(ref_index)
	path(known_sites_index)
	
	output:
	tuple val(id),
		path("${id}.recal.bam"),
		path("${id}.recal.bam.bai"),
		path("${id}.recal.before.table"),
		path("${id}.recal.after.table"),
		path("${id}.recalibration.csv")
	
	script:
	"""
	
	echo "Running BQSR for ${id}"
	
	# BEFORE
	gatk --java-options "-Xmx4g" BaseRecalibrator \
		-R ${params.reference} \
		-I ${marked_bam} \
		--known-sites ${params.known_sites} \
		-O ${id}.recal.before.table
		
	gatk --java-options "-Xmx4g" ApplyBQSR \
		-R ${params.reference} \
		-I ${marked_bam} \
		--bqsr-recal-file ${id}.recal.before.table \
		-O ${id}.recal.bam
	
	# AFTER	
	gatk --java-options "-Xmx4g" BaseRecalibrator \
		-R ${params.reference} \
		-I ${id}.recal.bam \
		--known-sites ${params.known_sites} \
		-O ${id}.recal.after.table
	
	gatk AnalyzeCovariates  \
		-before ${id}.recal.before.table \
		-after ${id}.recal.after.table \
		-csv ${id}.recalibration.csv

	gatk BuildBamIndex -I ${id}.recal.bam -O ${id}.recal.bam.bai
	"""
}

process haplotype_caller {
	
	tag "${id}_${interval_id}"
	
	container 'oras://community.wave.seqera.io/library/bcftools_bwa_fastp_gatk4_samtools:0fe683a5007bf05c'  
	
	input:
	tuple val(id), path(recal_bam), path(bai), val(interval_id), path(interval_file)
	path(ref_index)
	
	output:
	tuple val(id), val(interval_id),
		path("${id}.${interval_id}.g.vcf.gz"),
		path("${id}.${interval_id}.g.vcf.gz.tbi")
		
	script:
	"""
	echo "Running HaplotypeCaller for ${id} on interval ${interval_id}"
	
	gatk --java-options "-Xmx8g" HaplotypeCaller \
		-R ${params.reference} \
		-I ${recal_bam} \
		-L ${interval_file} \
		-O ${id}.${interval_id}.g.vcf.gz \
		-ERC GVCF \
		--native-pair-hmm-threads 4
	
	gatk IndexFeatureFile \
		-I ${id}.${interval_id}.g.vcf.gz
	
	echo "Finished HaplotypeCaller for ${id} on interval ${interval_id}"
	"""
}

process genomicsdb_import {

	tag "$interval_id"
	
	container 'oras://community.wave.seqera.io/library/bcftools_bwa_fastp_gatk4_samtools:0fe683a5007bf05c'  
	
	publishDir "${params.genomics_db}", mode: 'copy'
	
	input:
	tuple val(interval_id), path(gvcfs), path(tbis), path(interval_file)
	
	output:
	tuple val(interval_id),
		path("genomicsdb_${interval_id}")
	
	script:
	
	def workspace = "genomicsdb_${interval_id}"
	
    def workspace_flag = params.db_mode == 'update'
        ? "--genomicsdb-update-workspace-path ${workspace}"
        : "--genomicsdb-workspace-path ${workspace}"
		
	def vcf_args = gvcfs.collect { "-V ${it}" }.join(" ")
	
	"""
	echo "Running GenomicsDBImport for interval ${interval_id} (${params.db_mode})"
	
	gatk --java-options "-Xmx8g" GenomicsDBImport \
		${workspace_flag} \
		-L ${interval_file} \
		--interval-merging-rule OVERLAPPING_ONLY \
		${vcf_args}
		
	echo "Finished GenomicsDBImport for ${interval_id}"
    """
}

process genotype_gvcfs {

	tag "$interval_id"
	
	container 'oras://community.wave.seqera.io/library/bcftools_bwa_fastp_gatk4_samtools:0fe683a5007bf05c'  
	
	publishDir "${params.joint}", mode: 'copy'
	
	input:
	tuple val(interval_id), path(genomicsdb_workspace), path(interval_file)
	
	output:
	tuple val(interval_id),
		path("${params.vcf_prefix}.${interval_id}.vcf.gz"),
		path("${params.vcf_prefix}.${interval_id}.vcf.gz.tbi")
	
	script:
	
	"""
	echo "Genotyping interval ${interval_id}"
	
	gatk --java-options "-Xmx12g" GenotypeGVCFs \
		-R ${params.reference} \
		-V gendb://${genomicsdb_workspace} \
		-L ${interval_file} \
		-O ${params.vcf_prefix}.${interval_id}.vcf.gz
	
	echo "Done interval ${interval_id}"
	"""
}

process concat_vcfs {

	tag "final_cohort_merge"
	
	container 'oras://community.wave.seqera.io/library/bcftools_bwa_fastp_gatk4_samtools:0fe683a5007bf05c'  
	
	publishDir "${params.merged}", mode: 'copy'

	input:
	path(vcf_pairs)
	
	output:
	tuple path("${params.vcf_prefix}.merged.raw.vcf.gz"),
		  path("${params.vcf_prefix}.merged.raw.vcf.gz.tbi")	
	
	script:
	"""
	echo "Merging VCFs"
	
	ls *.vcf.gz | sort -V > vcfs.list
	
	bcftools concat -a -Oz \
		-o ${params.vcf_prefix}.merged.raw.vcf.gz \
		-f vcfs.list
	
	echo "Indexing final VCF" 
	
	bcftools index -t ${params.vcf_prefix}.merged.raw.vcf.gz
	
	"""
}






















