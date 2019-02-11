/*
 * Parent.cpp
 *
 *  Created on: 19 nov. 2018
 *      Author: Louison Fresnais, François Courtin
 *      Project: SMMB-ACO and Genetic Algorithm for epistasis detection
 *      Under the supervision of Christine Sinoquet(Nantes University)
 */

#include "Parent.h"
#include "Population.h"

using namespace std;

Parent::Parent(int nb_sol, int nb_parents, int len_pattern, double median, double P_selection,  double_matrix_type Mpop_geno) : nb_sol(nb_sol), nb_parents(nb_parents), len_pattern(len_pattern), median(median), P_selection(P_selection), Mpop_geno(Mpop_geno)
{
	int_matrix_type Mparents_temp (nb_parents, 1);
	Mparents = Mparents_temp;
}

Parent::~Parent() {
	// TODO Auto-generated destructor stub
}

/*
 * *************************************************************
 */

//Random selection of parents in the population
void Parent::parents_selection(){
	//cout << "Parents selection" << endl;
	vector<int> parents_list = {};
	int selected_parent;
	int compt = 0; //if more than 20 try to find a unique parent then get the parent even if not unique

	for (int i = 0; i < nb_parents; i++){
		//Pick a random number, if lower than the median of population g2 score, pick a "good" parent
		if ((rand()%100)+1 > P_selection){
			//Try to find a prent that is not already selected, 20 attempts or take the parent even if already selected
			do {
				selected_parent = rand()%nb_sol;
				compt ++;
				if (compt == 20){
					break;
				}
			}
			while (count(parents_list.begin(), parents_list.end(),selected_parent) != 0 or Mpop_geno(selected_parent, len_pattern) < median);
			parents_list.push_back(selected_parent);

		}else{
			//Pick a random number, if higher than the median of population g2 score, pick a "bad" parent
			//Try to find a prent that is not already selected, 20 attempts or take the parent even if already selected
			do {
				selected_parent = rand()%nb_sol;
				compt ++;
				if (compt == 20){
					break;
				}
			}
			while (count(parents_list.begin(), parents_list.end(),selected_parent) != 0 or Mpop_geno(selected_parent, len_pattern) >= median);

			parents_list.push_back(selected_parent);

		}
		Mparents(i, 0) = selected_parent;
	}
}

/*
 * *************************************************************
 */

int_matrix_type Parent::get_MParents(){
	return(Mparents);
}

/*
 * *************************************************************
 */

void Parent::display_parents(){
	cout << "parents sélectionnés" << Mparents << endl;
}
