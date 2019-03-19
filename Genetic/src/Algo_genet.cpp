/*
 * Algo_genet.cpp
 *
 *  Created on: 19 nov. 2018
 *      Author: Louison Fresnais, François Courtin
 *      Project: SMMB-ACO and Genetic Algorithm for epistasis detection
 *      Under the supervision of Christine Sinoquet(Nantes University)
 */

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/array.hpp>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <string>
#include "includes/Parametersfileparsing.hpp"
#include "includes/Population.h"
#include "includes/Parent.h"
#include "includes/ContingencyTable.h"
#include "includes/TheoricalTable.h"
#include "includes/G2test.h"
#include "includes/Child.h"
#include "includes/Output.h"
#include "includes/datainput.hpp"


using namespace std;


int main (int argc, char *argv[]) {
	srand (time(NULL));

	string parameter_file = argv[4];

	Parameters_file_parsing parameters(parameter_file);


	//Variable asignation from the parsing of the parameters files.
	char sep = parameters.sep;
	int nb_line_header = parameters.header_nrows;
	int len_pop = parameters.len_pop;
	int len_pattern = parameters.len_pattern;
	int nb_parents = parameters.nb_parents;
	double alpha = parameters.alpha;
	int nb_it = parameters.nb_it;
	double P_mutation = parameters.P_mutation;
	double P_selection = parameters.P_selection;
	string namefile = argv[1];
	string geno_path = argv[2];
	string pheno_path = argv[3];


	/*
	 * DATA READING AND POPULATION INITIALIZATION
	 */


	//Create a usefull genotype boost matrix from the input files
	Data_input datas_geno(geno_path, sep, nb_line_header);
	int_matrix_type Mgeno = datas_geno.read();

	//Create a usefull phenotype boost matrix from the input files
	Data_input datas_pheno(pheno_path, sep, nb_line_header);
	int_matrix_type Mpheno = datas_pheno.read();

	//Get the header from the input files into a vector
	Data_input header_line(geno_path, sep, nb_line_header);
	vector<string> header = header_line.get_snps();

	//Constructor for the population
	Population population(Mgeno, Mpheno, len_pop, len_pattern);

	//Check if the length of the population is not to high.
	population.check_pop_len();

	//Initialize the initial population of solutions
	population.init_pop_geno();
	//population.display_geno_sol();


	/*
	 * G2 TEST COMPUTED ON ALL THE INITIAL POPULATION
	 */


	int_matrix_type Msol_geno(1, len_pattern); //temp solution for the loop
	//Declaration of a vector for the G2 scores
 	vector<double> G2_res;

 	//Declaration and initialization of the non reliable G2 test count
 	int not_reliable_compt = 0;

 	//Loop to execute the G2 test for all the initial solutions, one solution at each loop turn
	for (int i = 0; i < len_pop; i++){
		//Get a solution from the population and put it in a temporary matrix
		for (int j = 0; j < len_pattern; j++){
			Msol_geno(0, j) = population.get_Mpop_geno()(i, j);
		}
		//Constructor for the contingency table
		ContingencyTable cont_table_pop(Mgeno, Mpheno, Msol_geno, len_pattern);
		//Set a list of all possible genotype patterns
		cont_table_pop.set_pattern_list();
		//Set the contingency table
		cont_table_pop.set_table();
		//cont_table_pop.display_table();

		//Constructor for the theorical contingency table
		TheoricalTable theo_table_pop(len_pattern, cont_table_pop.get_cont_table());
		//Set the theorical contingency table
		theo_table_pop.set_table();
		//theo_table_pop.display_table();

		//Constructor for the G2 test
		G2test G2_pop(cont_table_pop.get_cont_table(), theo_table_pop.get_theo_table());
		//Execute the G2 test on the choosen solution
		G2_pop.run_G2(not_reliable_compt, false);
		//G2_pop.display_g2();

		//Put the G2 score in the vector of the scores
		G2_res.push_back(G2_pop.get_g2());

		//Put the G2 score and the p-value of the solution into the population matrix according to the choosen solution
		population.set_Mpop_geno(i,len_pattern,G2_pop.get_g2());
		population.set_Mpop_geno(i,len_pattern+1,G2_pop.get_pval());

	}

	//Declare a variable for the median of the G2 scores
	double median;

	//Sort the G2 score vector
	sort(G2_res.begin(), G2_res.end());

	//Calculate the median of the G2 scores in the vector.
	median = G2_res[G2_res.size()/2];


	/*
	 * BEGINING OF THE LOOP OF N ITERATIONS
	 */


	int iterator = 0;
	//Begining of the n iterations of the methode
	while (iterator < nb_it){
		srand (time(NULL));

		/*
		 * PARENT SELECTION AND CROSSING OVER / MUTATION
		 */

		//Constructor for the parents
		Parent parents(len_pop, nb_parents, len_pattern, median, P_selection, population.get_Mpop_geno());
		//Select random paretns in the current population
		parents.parents_selection();
		//parents.display_parents();

		//Constructor for the children
		Child children(population.get_Mpop_geno(), parents.get_MParents(), len_pattern, P_mutation, Mgeno.size2());
		//Compute the crossing over between the parents to create the children
		children.set_children();
		//Compute the mutation (if there is one)
		children.mutation();
		//children.display_children();

		/*
		 * G2 TEST COMPUTED ON THE CHILDREN
		 */

		int_matrix_type Mchild(1, len_pattern); //Temp matrix for the loop

		//Loop to execute the G2 test for all the children made before
		for (int i = 0; i < nb_parents; i++){
			//Get a solution from the children and put it in a temporary matrix
			for (int j = 0; j < len_pattern; j++){
				Mchild(0, j) = children.get_MChildren()(i, j);
			}

			//Constructor for the contingency table for the choosen child
			ContingencyTable cont_table_child(Mgeno, Mpheno, Mchild, len_pattern);
			//Set all the possible genotype pattern
			cont_table_child.set_pattern_list();
			//Create the child contingency table
			cont_table_child.set_table();
			//cont_table_child.display_table();

			//Constructor for the Theorical contingency table for the child
			TheoricalTable theo_table_child(len_pattern, cont_table_child.get_cont_table());
			//Set the theorical contingency table of the child.
			theo_table_child.set_table();
			//theo_table_child.display_table();

			//Constructor for the G2 test
			G2test G2_child(cont_table_child.get_cont_table(), theo_table_child.get_theo_table());
			//Run the G2 test with the child contingency tables
			G2_child.run_G2(not_reliable_compt, true);
			//G2_child.display_g2();

			/*
			 * PARENTS SUBTITUTING BY BETTER CHILDREN
			 */

			//If the child has a better G2 score and lower p-value that its parent, replace the paretn by the child in the population matrix.
			if ((G2_child.get_g2() > population.get_Mpop_geno()(parents.get_MParents()(i,0), len_pattern)) and (G2_child.get_pval() < population.get_Mpop_geno()(parents.get_MParents()(i,0), len_pattern+1))){
				for (int j = 0; j < len_pattern; j++){
					population.set_Mpop_geno( parents.get_MParents()(i,0) , j , children.get_MChildren()(i,j) );
				}
				population.set_Mpop_geno( parents.get_MParents()(i,0) , len_pattern , G2_child.get_g2() );
				population.set_Mpop_geno( parents.get_MParents()(i,0) , len_pattern+1 , G2_child.get_pval() );
			}
		}


		//Clear the G2 score vector
		G2_res.clear();
		//set the new G2 score vector with the new scores
		for (int i = 0; i < len_pop; i++){
			G2_res.push_back(population.get_Mpop_geno()(i,len_pattern));
		}
		//Sort the G2 score vector
		sort(G2_res.begin(), G2_res.end());
		//Compute the new median
		median = G2_res[G2_res.size()/2]; //Median of the solutions' G2

		iterator++;
		//cout << iterator << "/";
	}

	/*
	 * END OF THE LOOP OF N ITERATIONS
	 */

	/*
	 * SOLUTIONS OUTPUT
	 */


	//Constructor for the output module
	Output output(population.get_Mpop_geno(), header, len_pattern, len_pop, namefile, alpha);

	//output.set_list_pattern();
	//Set the list of the final solutions sorting it by p-value
	output.set_best_sol();
	//Write the list in a file with the score and the p-value
	output.write_best_sol();

	//Display the number of non reliable G2 test
	//cout << "Number of not reliable G2 test : " << not_reliable_compt << endl;
	//cout << "END" << endl;
}
