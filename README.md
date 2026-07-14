# Advanced Crop Variant Discovery

A generalizable framework for advanced variant discovery in crop genomics.

## Status

Under active development(!)

This repository is currently being developed alongside the manuscript:

A Generalizable Framework for Advanced Variant Discovery in Crop Genomics.

Workflows, documentation, and software components may change until the first stable release.

## Authors

- Thomas Bergmann

## Citation

Citation information and DOI will be added upon publication.

## Description

A generalizable framework for advanced SNP and INDEL discovery in crop species that lack
comprehensive, high-confidence genomic reference resources.

The framework adapts GATK Best Practices for crop genomics by enabling Variant QUality Score Recalibration (VQSR)
without requiring predefined truth sets. Using unsupervised machine learning approaches, it generates reliable training
datasets for variant recalibration and filtering across diverse crop species and population sizes.

## Associated Manuscripts

This repository supports the manuscript:

**Bergmann T., MacNish T.R., Batley J., Edwards, D.**
* A Generalizable Framework for Avdanced Variant Discovery in Crop Genomics*

### Abstract

In crop genomics, variant discovery pipelines are constrained by the limited availability of standardised, high-confidence reference resources compared with human genomics. 
This complicates the adoption of GATK Best Practices in plant population studies despite their accuracy, flexibility, and scalability. Publicly available plant genomic resources are often based on outdated reference genomes, heterogenous pipelines, and generic hard-filtering thresholds derived from human studies, often not applicable to crop genomic data. 
Here, we present a crop-adapted variant recalibration framework that enables advanced SNP and insertion-deletion (INDEL) discovery without requiring predefined high-confidence reference variant sets. Using Brassica napus as a model system, we utilize the Illumina 60k Brassica SNP array to generate a reference-specific truth set for data pre-processing and to evaluate an unsupervised machine-learning approach for constructing variant quality score recalibration (VQSR) training datasets. Using the SNP array-derived truth set as a reference, we demonstrate that Gaussian mixture model (GMM)-based clustering effectively distinguishes high- and low-confidence variant profiles and provides a flexible approach for generating reliable training datasets for VQSR.
We further evaluate the framework in diploid Brassica oleracea and hexaploid wheat (Triticum aestivum), spanning populations that differ in size and genetic diversity. Across species and population scales, GMM-based clustering consistently generates effective VQSR training datasets and enables advanced filtering of both SNPs and INDELs without reliance on external reference resources. This adapted GATK Best Practices framework provides a significant advancement in variant discovery in crop species lacking comprehensive genomic resources.

## Key Features

- GATK-based variant discovery workflow
- Crop-adapted VQSR framework
- Unsupervised Gaussian Mixture Model (GMM) clustering
- SNP and INDEL recalibration without external truth resources
- Applicable across diverse crop species and ploidy levels
- Scalable to large population studies

## Species Evaluated

- *Brassica napus*
- *Brassica oleracea*
- *Triticum aestivum* (in progress)

## Citation

Citation information and DOI will be added upon publication.

The workflows implemented in this repository have been applied in"

**MacNish et al.**
*Single Nucleotide Polymorphisms and Haplotypes Provide Complementary Information for Genomic Selection in Brassica napus*
DOI: TBD

Feedback and suggestions are welcome.


