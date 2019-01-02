#!/bin/bash
#run with ./launch_genetic.sh   name_of_the_output_files   geno_path   pheno_path

#ex : ./launch_genetic.sh res Utils/data_simulation/test_perturb_all_files/test_perturbgenotypes1.txt Utils/data_simulation/test_perturb_all_files/test_perturbphenotypes1.txt 

for i in {1..20}
do
  ./Genetic/Debug/Genetic "$1_$i" "$2" "$3"
done
