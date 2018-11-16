#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/array.hpp>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>

using namespace boost::numeric::ublas;

typedef boost::numeric::ublas::matrix<int> int_matrix_type;
typedef boost::array<int_matrix_type, 20> Marray_type;
typedef boost::array<int_matrix_type, 4> MChild_array_type;


//Methode initiation population(matrix of genotypes, size of the population)
Marray_type init_pop(int_matrix_type Mgeno, int nb_sol){



	int nb_indiv = Mgeno.size1();
	int nb_snp = Mgeno.size2();

	int_matrix_type Mpop (nb_indiv, 3); //Matrix of a solution

	//Array of matrix = list of the solutions
	Marray_type sol_list;

	int nb_indiv_pop = Mpop.size1();
    //int nb_snp_pop = Mpop.size2();
    int indiv_select;


    for (int i = 0; i < nb_sol; i++){
    	for (int j = 0; j < nb_indiv_pop; j++){
        	indiv_select = rand()%((nb_indiv)); //Hasard selection of individus in the genotypes matrix
        	Mpop (j, 0) = Mgeno(indiv_select, rand()%(nb_snp)); //hasard selection of 3 snp from an individu to forme the solution
        	Mpop (j, 1) = Mgeno(indiv_select, rand()%(nb_snp));
        	Mpop (j, 2) = Mgeno(indiv_select, rand()%(nb_snp));
    	}
    	sol_list[i] = {Mpop}; //add the solution matrix in the list of solutions
    }

    return sol_list;

}


//Crossover between 2*2 solutions hasardly pick.
MChild_array_type cross_sol(Marray_type pop_sol){

	//Hasard selection of 2 pairs of prents
	int parent1 = rand()%(pop_sol.size());
	int parent2 = rand()%(pop_sol.size());
	int parent3 = rand()%(pop_sol.size());
	int parent4 = rand()%(pop_sol.size());

	//Control if parents are not the same
	while (parent2 == parent1){
		parent2 = rand()%(pop_sol.size());
	}

	while (parent3 == parent1 or parent3 == parent2){
		parent3 = rand()%(pop_sol.size());
	}

	while (parent4 == parent1 or parent4 == parent2 or parent4 == parent3){
		parent4 = rand()%(pop_sol.size());
	}

	//first pair of parents
	int_matrix_type MParent1(50, 2);
	MParent1 = pop_sol[parent1];
	std::cout << "parent 1 " << ":" << MParent1 << std::endl;
	int_matrix_type MParent2(50, 2);
	MParent2= pop_sol[parent2];
	std::cout << "parent 2 " << ":" << MParent2 << std::endl;

	//second pair of parents
	int_matrix_type MParent3(50, 2);
	MParent3= pop_sol[parent3];
	std::cout << "parent 3 " << ":" << MParent3 << std::endl;
	int_matrix_type MParent4(50, 2);
	MParent4 = pop_sol[parent4];
	std::cout << "parent 4 " << ":" << MParent4 << std::endl;

	//Initialization of cut places in parents
	int cut1 = rand()%(MParent1.size2()-1)+1;
	int cut2 = rand()%(MParent2.size2()-1)+1;
	int cut3 = rand()%(MParent3.size2()-1)+1;
	int cut4 = rand()%(MParent4.size2()-1)+1;
	std::cout << "cut1 " << ":" << cut1 << std::endl;
	std::cout << "cut2 " << ":" << cut2 << std::endl;
	std::cout << "cut3 " << ":" << cut3 << std::endl;
	std::cout << "cut4 " << ":" << cut4 << std::endl;

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

	std::vector<int> child1_part2;
	std::vector<int> child2_part2;
	std::vector<int> child3_part2;
	std::vector<int> child4_part2;

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


	std::cout << "CROSSING" << std::endl;
	std::cout << "child 1 :" << Mchild1 << std::endl;
	std::cout << "child 2 :" << Mchild2 << std::endl;
	std::cout << "child 3 :" << Mchild3 << std::endl;
	std::cout << "child 4 :" << Mchild4 << std::endl;


	//MUTATION GENERATION

	int mutation;

	if (double(rand()%101)/100 < 0.1){ //If proba < 0.1 then make a mutation
		std::cout << "MUTATION ON CHILD1" << std::endl;
		int j = rand()%int(Mchild1.size2()); //Set the snp to mutate
		for (int i = 0; i < int(Mchild1.size1()); i++){
			do{
				mutation = int(rand()%3); //Mutate every snp choosen in a different genotype
			}while (Mchild1(i, j) == mutation);
			Mchild1(i, j) = mutation;
		}
		std::cout << "child 1 :" << Mchild1 << std::endl;
	}

	if (double(rand()%101)/100 < 0.1){
		std::cout << "MUTATION ON CHILD2" << std::endl;
		int j = rand()%int(Mchild2.size2());
		for (int i = 0; i < int(Mchild2.size1()); i++){
			do{
				mutation = int(rand()%3);
			}while (Mchild2(i, j) == mutation);
			Mchild2(i, j) = mutation;
		}
		std::cout << "child 2 :" << Mchild2 << std::endl;
	}

	if (double(rand()%101)/100 < 0.1){
		std::cout << "MUTATION ON CHILD3" << std::endl;
		int j = rand()%int(Mchild3.size2());
		for (int i = 0; i < int(Mchild3.size1()); i++){
			do{
				mutation = int(rand()%3);
			}while (Mchild3(i, j) == mutation);
			Mchild3(i, j) = mutation;
		}
		std::cout << "child 3 :" << Mchild3 << std::endl;
	}

	if (double(rand()%101)/100 < 0.1){
		std::cout << "MUTATION ON CHILD4" << std::endl;
		int j = rand()%int(Mchild4.size2());
		for (int i = 0; i < int(Mchild4.size1()); i++){
			do{
				mutation = int(rand()%3);
			}while (Mchild4(i, j) == mutation);
			Mchild4(i, j) = mutation;
		}
		std::cout << "child 4 :" << Mchild4 << std::endl;
	}


	MChild_array_type Children_list;

	Children_list[0] = {Mchild1};
	Children_list[1] = {Mchild2};
	Children_list[2] = {Mchild3};
	Children_list[3] = {Mchild4};

	return(Children_list);

}

int main () {

	srand (time(NULL));
	int_matrix_type m (50, 10);

    int nb_indiv = m.size1();
    int nb_snp = m.size2();

    for (int i = 0; i < nb_indiv; ++ i)
        for (int j = 0; j < nb_snp; ++ j)
            m (i, j) = rand()%3;
    std::cout << "matrice de départ : " << m << std::endl;

    Marray_type Apop = init_pop(m, 20);

    for (int i = 0; i < 20; i++)
    	std::cout << "Solution " << i+1 << ":" << Apop[i] << std::endl;

    MChild_array_type AChildren = cross_sol(Apop);


}
