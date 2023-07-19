## This workflow is not intended to be directly excuted but to use it as a guide in terminal to do the analysis

#### For some software we recommend the use of conda

## To install conda:

## Download conda installer for Linux https://www.anaconda.com/download/ 

bash Anaconda-latest-Linux-x86_64.sh

## Despues crea los environments antes de instalar cualquier programa

conda create -n environment_name

conda activate environment_name

conda install -c bioconda environment_name


## Extract raw reads from a tar file

tar -xzvf strains.tar.gz


## Trimming the reads
## https://github.com/OpenGene/fastp

conda create -n fastp 

conda activate fastp 

conda install -c bioconda fastp

conda activate fastp

## Go to the directoy where raw reads have been extracted:

for R1 in *1.fq.gz
 do echo $R1
 R2=`echo $R1 | sed 's/_1.fq.gz/_2.fq.gz/'`
 NAME=`echo $R1 | sed 's/_E.\+//'`
 fastp --in1 $R1 --in2 $R2 --out1 "$NAME"_1_trimmed.fq.gz --out2 "$NAME"_2_trimmed.fq.gz -l 50 --dedup -w 8 --unpaired1 "$NAME"_unpaired.fq.gz --unpaired2 "$NAME"_unpaired.fq.gz 
done

conda deactivate ## important to go back to base with everything you have installed in /usr/bin or etc.

#Next step is to assemble genomes using unicycler v0.5 // No pilon polish
mkdir unicycler_output

for i in *_trimmed.fq.gz
 do 
 NAME1=`echo $i | grep "_1_trimmed"`
 NAME2=`echo $NAME1 | sed "s/_1_trimmed/_1_trimmed/"`
 NEWNAME=` echo $NAME1 | sed 's/__1_trimmed.\+//'`
 unicycler -1 $NAME1 -2 $NAME2 -o unicycler_output/$NEWNAME --verbosity 2 -t 20 --spades_options "--memory 16" 
done
 
## we end with assembled genomes in *.fasta