nextflow.enable.dsl=2


/*
 * =========================
 * PARAMETERS
 * =========================
 */

params.reference = "${workflow.projectDir}/reference/DarmorV10_Chromosomes_Only.fa"
params.raw_SNPs = "${workflow.projectDir}/results/joint_genotyping/merged/Samples_Bna_Minicore_1-135.merged.raw.SNPs.vcf.gz"
params.truth_set = "${workflow.projectDir}/results/joint_genotyping/merged/Samples_Bna_Minicore_1-135.TRUTH_SET.vcf.gz"
params.good_SNPs = "${workflow.projectDir}/training/VQSR_Training_Good.DV10.Minicore.vcf.gz"
params.medium_SNPs = "${workflow.projectDir}/training/VQSR_Training_Medium.DV10.Minicore.vcf.gz"
params.outdir = "${workflow.projectDir}/results/vqsr"
params.vcf_prefix   = "Samples_Bna_Minicore_1-135"

/*
 * =========================
 * PROCESS: VariantRecalibrator
 * =========================
 */


process VariantRecalibrator {
	
	tag "VQSR SNP"
	
	container 'oras://community.wave.seqera.io/library/bcftools_bwa_fastp_gatk4_samtools:0fe683a5007bf05c'
	
	publishDir "${ params.outdir }", mode: 'copy'
	
	input:

	path ref
	path ref_fai
	path ref_dict

	path raw_vcf
	path raw_index

	path truth
	path truth_index

	path good
	path good_index

	path medium
	path medium_index
	
	output:
	path "output.SNPs.minicore.recal"
	path "output.SNPs.minicore.tranches"
	path "output.SNPs.minicore.plots.R"
	
	script:
	"""
	
	gatk --java-options "-Xmx32g -Xms32g" VariantRecalibrator \
		-R ${ref} \
		-V ${raw_vcf} \
		--trust-all-polymorphic \
		-mode SNP \
		-an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum \
		--resource:GMM_good,known=false,training=true,truth=true,prior=12.0 ${good} \
		--resource:GMM_medium,known=false,training=true,truth=false,prior=10.0 ${medium} \
		--resource:truth_set,known=true,training=false,truth=false,prior=2.0 ${truth} \
		-O output.SNPs.minicore.recal \
		--tranches-file output.SNPs.minicore.tranches \
		-rscript-file output.SNPs.minicore.plots.R \
		--dont-run-rscript
	"""
}

/*
 * =========================
 * PROCESS: ApplyVQSR
 * =========================
 */

process ApplyVQSR {

	tag "Apply VQSR"
	
	publishDir "${ params.outdir }", mode: 'copy'

	container 'oras://community.wave.seqera.io/library/bcftools_bwa_fastp_gatk4_samtools:0fe683a5007bf05c'
	
	input:
	path ref
	path ref_fai
	path ref_dict

	path raw_vcf
	path raw_index

	path recal
	path tranches
	
	output:
	path "${params.vcf_prefix}.SNPs.VQSR.99.9.tranche.vcf.gz"
	
	script:
	"""
	
	gatk IndexFeatureFile \
		-I ${recal}

	gatk --java-options "-Xmx16g -Xms16g" ApplyVQSR \
		-R ${ref} \
		-V ${raw_vcf} \
		-O ${params.vcf_prefix}.SNPs.VQSR.99.9.tranche.vcf.gz \
		--truth-sensitivity-filter-level 99.9 \
		--recal-file ${recal} \
		--tranches-file ${tranches} \
		-mode SNP
	"""
}	

/*
 * =========================
 * WORKFLOW
 * =========================
 */

workflow {

	
	ref_ch        = Channel.fromPath(params.reference)
	ref_fai_ch    = Channel.fromPath(params.reference + ".fai")
	ref_dict_ch   = Channel.fromPath(params.reference.replace(".fa", ".dict"))

	raw_ch        = Channel.fromPath(params.raw_SNPs)
	raw_index_ch  = Channel.fromPath(params.raw_SNPs + ".tbi")
	
	truth_ch      = Channel.fromPath(params.truth_set)
	truth_idx_ch  = Channel.fromPath(params.truth_set + ".tbi")

	good_ch       = Channel.fromPath(params.good_SNPs)
	good_idx_ch   = Channel.fromPath(params.good_SNPs + ".tbi")

	medium_ch     = Channel.fromPath(params.medium_SNPs)	
	medium_idx_ch = Channel.fromPath(params.medium_SNPs + ".tbi")


	recal_out = VariantRecalibrator(
		ref_ch,
		ref_fai_ch,
		ref_dict_ch,
		raw_ch,
		raw_index_ch,
		truth_ch,
		truth_idx_ch,
		good_ch,
		good_idx_ch,
		medium_ch,
		medium_idx_ch
	)
	
	ApplyVQSR(
		ref_ch,
		ref_fai_ch,
		ref_dict_ch,
		raw_ch,
		raw_index_ch,
		recal_out[0],
		recal_out[1]
	)
}
	
	
	
	
	
	
	
	
