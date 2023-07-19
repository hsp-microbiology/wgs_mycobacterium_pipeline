############################# Identification using assembled fastas ###################################

## pubMLST Species using python script copied from https://pubmlst.org/species-id/species-identification-via-api 
for i in BCN*.fasta
 do
 NAME=`echo $i | sed 's/.fasta//g'`
 python ~/pubmlst_species_identification.py -f $i > "$NAME".report
done

## Summarize reports into csv

for i in *.report
 do 
 NAME=`echo $i | sed 's/.report//g'`
 echo "name:$NAME" >> $i
 done

cat *.report > all_reports


 ## OrthoANI https://www.ezbiocloud.net/tools/orthoani 

java -jar ~/path_to_OUA/OAU.jar -u ~/path_to_OUA/usearch --fastadir ~/working_directory/ --out ANI_genome_comparison

