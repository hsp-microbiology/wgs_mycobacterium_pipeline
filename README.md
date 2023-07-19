# wgs_mycobacterium_pipeline
This is a repository to upload the bioinformatics pipeline used for the Whole Genome Sequencing analysis of Mycobacterium avium complex isolates. This repository contains scripts in bash and R 
for processing raw Illumina NGS data. This is the workflow I follow for Species Identification and Genomic epidemiology using available programs for Linux. Briefly, we begin with raw sequencing data,
then trimmed reads, to quality control. Then with trimmed reads of good quality you can perform downstream analysis.
## De Novo assembly
Using Unicycler, and then downstream analysis using prokka, Roary and RAxML.
## Map to Reference
Using Snippy, Gubbins, snp-sites, RAxML and R.
## Direct Cluster Analysis
Using SKA

## Disclosure
This code works well for our studies but is subjected to changes in the future. These scripts show our workflow and are just wrappers that use different softwares available. These scripts are not 
meant to be executed but rather to use as guides when working on terminal.
