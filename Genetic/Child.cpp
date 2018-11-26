/*
 * Child.cpp
 *
 *  Created on: 19 nov. 2018
 *      Author: courtin
 */
#include "Child.h"

Child::Child(float_matrix_type Mpop_geno, int_matrix_type MParents, int len_pattern) : Mpop_geno(Mpop_geno), MParents(MParents), len_pattern(len_pattern)
{
	int_matrix_type MChlidren_temp(MParents.size1(), len_pattern);
	MChildren = MChlidren_temp;
}

Child::~Child() {
	// TODO Auto-generated destructor stub
}

void Child::set_children(){
	int it = 0;


	if (len_pattern < 3){
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

void Child::mutation(){
	for (int i = 0; i < MChildren.size1(); i++){
		if ((rand()%100)+1 < 10){
			MChildren(i, rand()%len_pattern) = rand()%Mpop_geno.size1();
			cout << "mutation on child " << i << endl;
		}
	}
}

void Child::display_children(){
	cout << "children" << MChildren << endl;
}

int_matrix_type Child::get_MChildren(){
	return(MChildren);
}
