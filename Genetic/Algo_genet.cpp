#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/array.hpp>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>

using namespace std;

typedef boost::numeric::ublas::matrix<int> int_matrix_type;
typedef boost::numeric::ublas::matrix<int_matrix_type> matrix_of_int_matrix_type;


//Methode initiation population(matrix of genomatrice de départ types, size of the population)
matrix_of_int_matrix_type init_pop(int_matrix_type Mgeno, int_matrix_type Mpheno, int nb_sol){


	int nb_indiv = Mgeno.size1();
	int nb_snp = Mgeno.size2();

	int_matrix_type Mpop (nb_indiv, 3); //Matrix of a solution
	int_matrix_type Mpop_pheno(nb_indiv, 1);

	//Array of matrix = list of the solutions
	matrix_of_int_matrix_type sol_list(2, nb_sol);

	int nb_indiv_pop = Mpop.size1();
    //int nb_snp_pop = Mpop.size2();
    int indiv_select;


    for (int i = 0; i < nb_sol; i++){
    	for (int j = 0; j < nb_indiv_pop; j++){
        	indiv_select = rand()%((nb_indiv)); //Hasard selection of individus in the genotypes matrix
        	Mpop (j, 0) = Mgeno(indiv_select, rand()%(nb_snp)); //hasard selection of 3 snp from an individu to forme the solution
        	Mpop (j, 1) = Mgeno(indiv_select, rand()%(nb_snp));
        	Mpop (j, 2) = Mgeno(indiv_select, rand()%(nb_snp));
        	Mpop_pheno (j, 0) = Mpheno(indiv_select, 0);
    	}
    	sol_list(0, i) = Mpop; //add the solution matrix in the list of solutions
    	sol_list(1, i) = Mpop_pheno;
    }

    return sol_list;

}


int_matrix_type parent_selection(matrix_of_int_matrix_type Mpop){


	//Hasard selection of 2 pairs of prents
	int parent1 = rand()%(Mpop.size2());
	int parent2 = rand()%(Mpop.size2());
	int parent3 = rand()%(Mpop.size2());
	int parent4 = rand()%(Mpop.size2());

	//Control if parents are not the same
	while (parent2 == parent1){
		parent2 = rand()%(Mpop.size2());
	}

	while (parent3 == parent1 or parent3 == parent2){
		parent3 = rand()%(Mpop.size2());
	}

	while (parent4 == parent1 or parent4 == parent2 or parent4 == parent3){
		parent4 = rand()%(Mpop.size2());
	}

	int_matrix_type Mparents(1, 4);
	Mparents(0, 0) = parent1;
	Mparents(0, 1) = parent2;
	Mparents(0, 2) = parent3;
	Mparents(0, 3) = parent4;

	return(Mparents);

}

//Crossover between 2*2 solutions hasardly pick.
matrix_of_int_matrix_type cross_sol(matrix_of_int_matrix_type Mpop, int_matrix_type Mparents){


	//first pair of parents
	int_matrix_type MParent1(50, 2);
	MParent1 = Mpop(0, Mparents(0, 0));
	cout << "parent 1 " << ":" << MParent1 << endl;
	int_matrix_type MParent2(50, 2);
	MParent2= Mpop(0, Mparents(0, 1));
	cout << "parent 2 " << ":" << MParent2 << endl;

	//second pair of parents
	int_matrix_type MParent3(50, 2);
	MParent3= Mpop(0, Mparents(0, 2));
	cout << "parent 3 " << ":" << MParent3 << endl;
	int_matrix_type MParent4(50, 2);
	MParent4 = Mpop(0, Mparents(0, 3));
	cout << "parent 4 " << ":" << MParent4 << endl;

	//Initialization of cut places in parents
	int cut1 = rand()%(MParent1.size2()-1)+1;
	int cut2 = rand()%(MParent2.size2()-1)+1;
	int cut3 = rand()%(MParent3.size2()-1)+1;
	int cut4 = rand()%(MParent4.size2()-1)+1;
	cout << "cut1 " << ":" << cut1 << endl;
	cout << "cut2 " << ":" << cut2 << endl;
	cout << "cut3 " << ":" << cut3 << endl;
	cout << "cut4 " << ":" << cut4 << endl;

	int_matrix_type Mchild1(MParent1.size1(), cut1 + MParent2.size2() - cut2);
	int_matrix_type Mchild2(MParent2.size1(), cut2 + MParent1.size2() - cut1);
	int_matrix_type Mchild3(MParent3.size1(), cut3 + MParent4.size2() - cut4);
	int_matrix_type Mchild4(MParent4.size1(), cut4 + MParent3.size2() - cut3);


	//Settinf first part of the children
	for (int i = 0; i < int(Mchild1.size1()); i++){
		for (int j = 0; j < cut1; j++){
			Mchild1(i, j) = MParent1(i, j);
		}
	}

	for (int i = 0; i < int(Mchild2.size1()); i++){
		for (int j = 0; j < cut2; j++){
			Mchild2(i, j) = MParent2(i, j);
		}
	}

	for (int i = 0; i < int(Mchild3.size1()); i++){
		for (int j = 0; j < cut3; j++){
			Mchild3(i, j) = MParent3(i, j);
		}
	}

	for (int i = 0; i < int(Mchild4.size1()); i++){
		for (int j = 0; j < cut4; j++){
			Mchild4(i, j) = MParent4(i, j);
		}
	}


	//Setting seconde part of the children in vectors
	//adding vectors on the children to complete the crossover

	vector<int> child1_part2;
	vector<int> child2_part2;
	vector<int> child3_part2;
	vector<int> child4_part2;

	int vect_it;

	for (int i = 0; i < int(MParent2.size1()); i++){
		child1_part2.clear();
		vect_it = 0;
		for (int j = cut2; j < int(MParent2.size2()); j++){ //Putting second part on a vector
			child1_part2.push_back(MParent2(i, j));
		}
		for (int j = cut1; j < int(Mchild1.size2()); j++){ //Adding the vector to the chil matrix
			Mchild1(i, j) = child1_part2.at(vect_it);
			vect_it ++;
		}
	}

	for (int i = 0; i < int(MParent1.size1()); i++){
		child2_part2.clear();
		vect_it = 0;
		for (int j = cut1; j < int(MParent1.size2()); j++){
			child2_part2.push_back(MParent1(i, j));
		}
		for (int j = cut2; j < int(Mchild2.size2()); j++){
			Mchild2(i, j) = child2_part2.at(vect_it);
			vect_it ++;
		}
	}

	for (int i = 0; i < int(MParent4.size1()); i++){
		child3_part2.clear();
		vect_it = 0;
		for (int j = cut4; j < int(MParent4.size2()); j++){
			child3_part2.push_back(MParent4(i, j));
		}
		for (int j = cut3; j < int(Mchild3.size2()); j++){
			Mchild3(i, j) = child3_part2.at(vect_it);
			vect_it ++;
		}
	}

	for (int i = 0; i < int(MParent3.size1()); i++){
		child4_part2.clear();
		vect_it = 0;
		for (int j = cut3; j < int(MParent3.size2()); j++){
			child4_part2.push_back(MParent3(i, j));
		}
		for (int j = cut4; j < int(Mchild4.size2()); j++){
			Mchild4(i, j) = child4_part2.at(vect_it);
			vect_it ++;
		}
	}


	cout << "CROSSING" << endl;
	cout << "child 1 :" << Mchild1 << endl;
	cout << "child 2 :" << Mchild2 << endl;
	cout << "child 3 :" << Mchild3 << endl;
	cout << "child 4 :" << Mchild4 << endl;


	//MUTATION GENERATION

	int mutation;

	if (double(rand()%101)/100 < 0.1){ //If proba < 0.1 then make a mutation
		cout << "MUTATION ON CHILD1" << endl;
		int j = rand()%int(Mchild1.size2()); //Set the snp to mutate
		for (int i = 0; i < int(Mchild1.size1()); i++){
			do{
				mutation = int(rand()%3); //Mutate every snp choosen in a different genotype
			}while (Mchild1(i, j) == mutation);
			Mchild1(i, j) = mutation;
		}
		cout << "child 1 :" << Mchild1 << endl;
	}

	if (double(rand()%101)/100 < 0.1){
		cout << "MUTATION ON CHILD2" << endl;
		int j = rand()%int(Mchild2.size2());
		for (int i = 0; i < int(Mchild2.size1()); i++){
			do{
				mutation = int(rand()%3);
			}while (Mchild2(i, j) == mutation);
			Mchild2(i, j) = mutation;
		}
		cout << "child 2 :" << Mchild2 << endl;
	}

	if (double(rand()%101)/100 < 0.1){
		cout << "MUTATION ON CHILD3" << endl;
		int j = rand()%int(Mchild3.size2());
		for (int i = 0; i < int(Mchild3.size1()); i++){
			do{
				mutation = int(rand()%3);
			}while (Mchild3(i, j) == mutation);
			Mchild3(i, j) = mutation;
		}
		cout << "child 3 :" << Mchild3 << endl;
	}

	if (double(rand()%101)/100 < 0.1){
		cout << "MUTATION ON CHILD4" << endl;
		int j = rand()%int(Mchild4.size2());
		for (int i = 0; i < int(Mchild4.size1()); i++){
			do{
				mutation = int(rand()%3);
			}while (Mchild4(i, j) == mutation);
			Mchild4(i, j) = mutation;
		}
		cout << "child 4 :" << Mchild4 << endl;
	}


	matrix_of_int_matrix_type Children_list(1, 4);

	Children_list(0, 0) = Mchild1;
	Children_list(0, 1) = Mchild2;
	Children_list(0, 2) = Mchild3;
	Children_list(0, 3) = Mchild4;

	return(Children_list);

}



//Fonction will compare children with parents and update the population
void contingency_table(int_matrix_type MChildren_pheno,int_matrix_type MChildren_geno){


	int_matrix_type table(2, 3);
	for (int i = 0; i < int(table.size1()); i++){
		for (int j = 0; i < int(table.size2()); i++){
			table(i, j) = 0;
		}
	}

	for (int i = 0; i < int(MChildren_pheno.size1()); i++){
		for (int j = 0; j < int(MChildren_geno.size2()); j++){
			if (MChildren_pheno(i, 0) == 0){
				if(MChildren_geno(i, j) == 0){
					table(0, 0) += 1;
				}
				if(MChildren_geno(i, j) == 1){
					table(0, 1) += 1;
				}
				if(MChildren_geno(i, j) == 2){
					table(0, 2) += 1;
				}
			}else{
				if(MChildren_geno(i, j) == 0){
					table(1, 0) += 1;
				}
				if(MChildren_geno(i, j) == 1){
					table(1, 1) += 1;
				}
				if(MChildren_geno(i, j) == 2){
					table(1, 2) += 1;
				}
			}
		}
	}
	cout << "table contingence : " << table << endl;
	//return(table);
}



int main () {

	srand (time(NULL));
	int_matrix_type Mgeno (50, 10);
	int_matrix_type Mpheno (50, 1);

    int nb_indiv_geno = Mgeno.size1();
    int nb_snp_pheno = Mgeno.size2();

    for (int i = 0; i < nb_indiv_geno; ++ i)
        for (int j = 0; j < nb_snp_pheno; ++ j)
        	Mgeno (i, j) = rand()%3;
    cout << "matrice de départ génotypes : " << Mgeno << endl;

    int nb_indiv_pheno = Mpheno.size1();

    for (int i = 0; i < nb_indiv_pheno; ++ i)
    	Mpheno (i, 0) = rand()%2;
    cout << "matrice de départ phénotypes : " << Mpheno << endl;

    int nb_indiv_pop = 20;

    matrix_of_int_matrix_type Mpop_pheno_geno(2, nb_indiv_pop);
    Mpop_pheno_geno = init_pop(Mgeno, Mpheno, nb_indiv_pop);

    matrix_of_int_matrix_type Mpop_geno(1, nb_indiv_pop);
    int_matrix_type Mpop_pheno(1, nb_indiv_pop);

    for (int i = 0; i < int(Mpop_pheno_geno.size2()); i++){
    	Mpop_geno(0, i) = Mpop_pheno_geno(0, i);
    }

    for (int i = 0; i < int(Mpop_pheno_geno.size2()); i++){
        	Mpop_pheno(0, i) = Mpop_pheno_geno(1, i);
        }

    for (int i = 0; i < nb_indiv_pop; i++)
    	cout << " genno solution " << i+1 << ":" << Mpop_geno(0, i) << endl;
    for (int i = 0; i < nb_indiv_pop; i++)
    	cout << "pheno solution " << i+1 << ":" << Mpop_pheno(0, i) << endl;


    int_matrix_type Mparents(1, 4);
    Mparents = parent_selection(Mpop_geno);

    matrix_of_int_matrix_type MChildren = cross_sol(Mpop_geno, Mparents);
    int_matrix_type cont_table(2, 3);
    int_matrix_type test(1, 1);
    test = MChildren(0, 0);

    contingency_table(Mpop_pheno, test);

}
//Ma matrice de phénotypes est une matrice de matrice il faut que je la change en matrice d'entiers.



