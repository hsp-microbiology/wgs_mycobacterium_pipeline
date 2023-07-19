####### map2ref ######

## avium ## We use ref H87 CP018363 from genbank
## chimaera ## We use ref CP015278.1 from genbank
## intracellulare ## We use BCN_8 assembled genome from our study

## To perform the map2ref and generate the fasta alignment we use snippy from https://github.com/tseemann/snippy 
## Here we present the analysis for M. avium, for chimaera or intracellulare you need to change the reference.fasta for each reference genome
sh runme.sh

snippy-clean_full_aln core.full.aln > clean_full.aln

## Remove recombination events usign Gubbins https://github.com/nickjcroucher/gubbins 
conda activate gubbins

run_gubbins.py --verbose --filter-percentage 99 -p gubbins --threads 20 clean_full.aln

conda deactivate 

## Call only SNPs using snp-sites https://github.com/sanger-pathogens/snp-sites 
snp-sites -c gubbins.filtered_polymorphic_sites.fasta > clean.core.aln

## Then we can infer the phylogenetic tree using Maximum-Likelyhood with RAxML https://cme.h-its.org/exelixis/web/software/raxml/ 
raxmlHPC-PTHREADS-AVX -s clean.core.aln -n mavium -T 20 -m GTRCAT -p 5000 -x 5000 -N 1000 -f a

## To visualize and annotate the tree we like to use figtree http://tree.bio.ed.ac.uk/software/figtree/ 
java -jar ~/path_to_executable/figtree.jar
## You can export pdf files from figtree or change the tree format for downstream analysis


############################# Coverage and other statistics #################################

## get coverage 
for i in BCN_*
 do 
 echo $i
 samtools mpileup $i/*.bam | awk 'BEGIN{C=0}; {C=C+$4}; END{print C "\t" C/NR}' > coverage/"$i"_total_and_average_coverage
 samtools mpileup $i/*.bam | awk -v X="${MIN_COVERAGE_DEPTH}" '$9>=X' | wc -l > coverage/"$i"_coverage_depth
 done
 
 ## Get statistics for alignment
 
for i in BCN_*
 do 
 echo $i
 samtools flagstats $i/*.bam > coverage/"$i"_stats
done

## Here we go from trimmed reads to a phylogenetic tree inferred from a map to reference alignment
