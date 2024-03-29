/*
 * Child.cpp
 *
 *  Created on: 19 nov. 2018
 *      Author: Louison Fresnais, François Courtin
 *      Project: SMMB-ACO and Genetic Algorithm for epistasis detection
 *      Under the supervision of Christine Sinoquet(Nantes University)
 */
#include "includes/Child.h"

Child::Child(double_matrix_type Mpop_geno, int_matrix_type MParents, int len_pattern, double P_mutation, int nb_snp) : Mpop_geno(Mpop_geno), MParents(MParents), len_pattern(len_pattern), P_mutation(P_mutation), nb_snp(nb_snp)
{
	int_matrix_type MChlidren_temp(MParents.size1(), len_pattern);
	MChildren = MChlidren_temp;
}

Child::~Child() {
	// TODO Auto-generated destructor stub
}

/*
 * *************************************************************
 */

//Execute the crossing over between the parents
void Child::set_children(){
	//cout << "Crossing over..." << endl;
	int it = 0;


	if (len_pattern == 2){
		while (it < int(MParents.size1())){
			MChildren(it,0) = Mpop_geno(MParents(it,0), 0);
			MChildren(it,1) = Mpop_geno(MParents(it+1,0),1);
			MChildren(it+1,0) = Mpop_geno(MParents(it+1,0),0);
			MChildren(it+1,1) = Mpop_geno(MParents(it,0),1);
			it = it + 2;
		}
	}else if(len_pattern == 3){
		while (it < int(MParents.size1())){
			MChildren(it,0) = Mpop_geno(MParents(it,0), 0);
			MChildren(it,1) = Mpop_geno(MParents(it+1,0),1);
			MChildren(it,2) = Mpop_geno(MParents(it+1,0),2);
			MChildren(it+1,0) = Mpop_geno(MParents(it+1,0),0);
			MChildren(it+1,1) = Mpop_geno(MParents(it,0),1);
			MChildren(it+1,2) = Mpop_geno(MParents(it,0),2);
			it = it + 2;
		}
	}

}

/*
 * *************************************************************
 */

//Execute the mutation on a child
void Child::mutation(){
	//cout << "Mutation..." << endl;
	vector<int> muted_child; //To sort the new child
	//Loop through all the created children
	for (int i = 0; i < int(MChildren.size1()); i++){
		//Pick a random number, if it's lower than the mutation probabilitu, compute the mutation
		if ((rand()%100)+1 < P_mutation){
			int muted_snp;
			int new_snp;
			do{
				muted_snp = rand()%len_pattern; //Select the snp to be mutated in the solution
				new_snp = rand()%nb_snp; //Select the new snp from all the possibilities
			}while (MChildren(i, muted_snp) == new_snp); //Check if the new snp isn't the same as the old one
			MChildren(i, muted_snp) = new_snp;
			for (int j = 0; j < len_pattern; j++){
				muted_child.push_back(MChildren(i, j));
			}
			//Sort the mutated child
			sort(muted_child.begin(), muted_child.end());
			//Put the new child into the child matrix
			for (int j = 0; j < len_pattern; j++){
				MChildren(i, j) = muted_child[j];
			}
		}
	}
}

/*
 * *************************************************************
 */

void Child::display_children(){
	//cout << "children" << MChildren << endl;
}

/*
 * *************************************************************
 */

int_matrix_type Child::get_MChildren(){
	return(MChildren);
}
