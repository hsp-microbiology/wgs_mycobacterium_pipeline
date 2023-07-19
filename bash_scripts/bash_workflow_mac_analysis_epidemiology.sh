################ EPIDEMIOLOGY ###########################

## Calculate the joint ancestral state reconstruction algorithm of Pupko et al using pyjar https://github.com/simonrharris/pyjar 

pyjar.py -a clean_core.aln -t mycobacterium_SNP.tre -o pyjar -v ## In order to use rPinecone or fastbaps in R

#### Create Sequence cluster by split kmer analysis using SKA https://github.com/simonrharris/SKA 

conda activate ska

for R1 in *_1_trimmed.fq.gz
 do 
 echo $i
 R2=`echo $R1 | sed 's/_1_/_2_/'`
 NAME=`echo $R1 | sed 's/_1_trimmed.fq.gz//'`
 ska fastq -o $NAME $R1 $R2
 done
 
ska summary *.skf

## Merge kmers and calculate pairwise distances

ska merge -o merged *.skf

ska summary all_samples.skf > summary

