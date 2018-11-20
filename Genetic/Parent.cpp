/*
 * Parent.cpp
 *
 *  Created on: 19 nov. 2018
 *      Author: courtin
 */

#include "Parent.h"
#include "Population.h"
#include <set>

Parent::Parent(int_matrix_type Mgeno, int_matrix_type Mpheno, int len_pop, int len_pattern, int nb_parents) : Population(Mgeno, Mpheno ,len_pop, len_pattern), nb_parents(nb_parents)
{
	int_matrix_type Mparents_temp (nb_parents, 1);
	Mparents = Mparents_temp;
}

Parent::~Parent() {
	// TODO Auto-generated destructor stub
}

void Parent::parents_selection(){
	//Hasard selection of n pairs of parents
	set<int> parents_list = {999};
	int selected_parent;
	for (int i = 0; i < nb_parents; i++){
		selected_parent = 999;
		while (parents_list.count(selected_parent) != 0){
			selected_parent = rand()%nb_sol;
		}
		parents_list.emplace(selected_parent);
		Mparents(i, 0) = selected_parent;
	}
}

void Parent::display_parents(){
	cout << "parents sélectionnés" << Mparents << endl;
}
