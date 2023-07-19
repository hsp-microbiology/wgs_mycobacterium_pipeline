## prokka annotation from unicycler output
## prokka https://github.com/tseemann/prokka 
## we created a custom database for mycobacteria
prokka-genbank_to_fasta_db Mycobacterium1.gdb Mycobacterium2.gdb ... > mycobacterium.faa

cd-hit -i mycobacterium.faa -o mycobacterium -T 0 -M 0 -g 1 -s 0.8 -c 0.9

rm -fv mycobacterium.faa mycobacterium.bak.clstr mycobacterium.clstr

makeblastdb -dbtype prot -in mycobacterium

mv mycobacterium.p* ~/path_to_prokka/prokka/db/genus/

## run prokka
for i in *.fasta
 do 
 NAME=`echo $i | sed 's/.fasta//'`
 prokka --outdir ~/prokka_output/$NAME --prefix $NAME --genus Mycobacterium --protein ~/path_to_prokka/prokka/db/genus/mycobacterium.faa --strain $NAME --cpus 0 "$i"
done

### Put all the .gff files in the same directory and run roary
## Roary: Pangenome analysis. From the *.gff files created from prokka

for i in BCN*/
 do echo $i
 cp "$i"*.gff ~/my_working_directory
done
 
# Check roary
roary -a

# Run Roary
roary -p 16 -v -e -n -i 90 -f roary *.gff

## Infer tree
cd roary 

raxmlHPC -s core_gene_alignment.aln -n roary_ntm_all_DATE -T 16 -m GTRCAT -p 5000 -x 5000 -N 99 -f a

## visualize tree
java -jar ~/path_to_figtree/figtree.jar
