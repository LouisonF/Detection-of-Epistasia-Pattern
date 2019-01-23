#!/bin/bash
#run with ./launch_genetic.sh   name_of_the_output_files   geno_path   pheno_path

#ex : ./launch_genetic.sh res /home/courtin/Documents/M2/ProjetC/Simu_naive_2snp_0.25 /home/courtin/Documents/M2/ProjetC/Simu_naive_2snp_0.25_pheno

_debut=$(date +%s)

for j in {1..100}
do
  for i in {1..20}
  do
    debut=$(date +%s)
    echo "FILE NUMBER : $j"
    echo "ITERATION NUMBER : $j"
    ./Genetic/Debug/Genetic "$1_$i" "$2/Naif_$j""_Genotype.txt" "$3/Naif_$j""_Phenotype.txt"
    fin=$(date +%s)

    duree=$(( $fin - $debut ))
    echo "$duree" >> tps.txt
  done
done
_fin=$(date +%s)

duree=$(( $_fin - $_debut ))
echo "$duree" >> tps.txt
