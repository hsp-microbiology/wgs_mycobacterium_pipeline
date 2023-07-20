# wgs_mycobacterium_pipeline
This is a repository to upload the bioinformatics pipeline used for the Whole Genome Sequencing analysis of Mycobacterium avium complex isolates. This repository contains scripts in bash and R 
for processing raw Illumina NGS data. This is the workflow I follow for Species Identification and Genomic epidemiology using available programs for Linux. Briefly, we begin with raw sequencing data,then trimmed reads, to quality control. Then with trimmed reads of good quality you can perform downstream analysis.
## De Novo assembly
Using Unicycler, and then downstream analysis using prokka, Roary and RAxML.
## Map to Reference
Using Snippy, Gubbins, snp-sites, RAxML and R.
## Direct Cluster Analysis
Using SKA
## Disclosure
This code works well for our studies but is subjected to changes in the future. These scripts show our workflow and are just wrappers that use different softwares available. These scripts are not 
meant to be executed but rather to use as guides when working on terminal.

### Software used in this pipeline
Programs:
fastp v0.23.2 (available at  https://github.com/OpenGene/fastp). 
Unicycler v0.5.0 (available at https://github.com/rrwick/Unicycler). 
SPAdes v3.15.3 (available at https://cab.spbu.ru/software/spades/). 
prokka v1.14.6 (available at https://github.com/tseemann/prokka). 
Roary v1.7.8 (available at https://github.com/sanger-pathogens/Roary). 
RAxML v8.2.12 (available at https://cme.h-its.org/exelixis/web/software/raxml/). 
figtree v1.4.4 (available at https://github.com/rambaut/figtree). 
Snippy v4.6.0 (available at  https://github.com/tseemann/snippy). 
Freebayes v1.3.2 (available at https://github.com/freebayes/freebayes).
Gubbins v3.2.1 (available at https://github.com/nickjcroucher/gubbins). 
snp-sites v2.5.1 (available at  https://github.com/sanger-pathogens/snp-sites).  
R v4.3.1 (available at https://cran.r-project.org/). 
rPinecone (available at  https://github.com/alexwailan/rpinecone).
SKA (available at  https://github.com/simonrharris/SKA).
Online resources:
PubMLST  (available at https://pubmlst.org/bigsdb?db=pubmlst_rmlst_seqdef_kiosk).
OrthoANI (available at https://www.ezbiocloud.net/tools/orthoani).
