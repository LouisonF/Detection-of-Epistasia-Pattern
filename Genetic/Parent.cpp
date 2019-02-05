/*
 * Parent.cpp
 *
 *  Created on: 19 nov. 2018
 *      Author: courtin
 */

#include "Parent.h"
#include "Population.h"


Parent::Parent(int nb_sol, int nb_parents, int len_pattern, double median, double P_selection,  double_matrix_type Mpop_geno) : nb_sol(nb_sol), nb_parents(nb_parents), len_pattern(len_pattern), median(median), P_selection(P_selection), Mpop_geno(Mpop_geno)
{
	int_matrix_type Mparents_temp (nb_parents, 1);
	Mparents = Mparents_temp;
}

Parent::~Parent() {
	// TODO Auto-generated destructor stub
}

void Parent::parents_selection(){
	//Hasard selection of n pairs of parents
	cout << "Parents selection" << endl;
	vector<int> parents_list = {};
	int selected_parent;
	for (int i = 0; i < nb_parents; i++){
		if ((rand()%100)+1 > P_selection){
			do {
				selected_parent = rand()%nb_sol;
			}
			while (count(parents_list.begin(), parents_list.end(),selected_parent) != 0 or Mpop_geno(selected_parent, len_pattern) < median);
			parents_list.push_back(selected_parent);
		}else{
			do {
				selected_parent = rand()%nb_sol;
			}
			while (count(parents_list.begin(), parents_list.end(),selected_parent) != 0 or Mpop_geno(selected_parent, len_pattern) >= median);
			parents_list.push_back(selected_parent);
		}
		Mparents(i, 0) = selected_parent;
	}
}

int_matrix_type Parent::get_MParents(){
	return(Mparents);
}

void Parent::display_parents(){
	cout << "parents sélectionnés" << Mparents << endl;
}
