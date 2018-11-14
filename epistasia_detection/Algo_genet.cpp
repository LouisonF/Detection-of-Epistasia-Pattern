#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/array.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

using namespace boost::numeric::ublas;

typedef boost::numeric::ublas::matrix<int> int_matrix_type;
typedef boost::array<int_matrix_type, 20> Marray_type;


//Methode initiation population(matrix of genotypes, size of the population)
Marray_type init_pop(int_matrix_type Mgeno, int nb_sol){



	int nb_indiv = Mgeno.size1();
	int nb_snp = Mgeno.size2();

	int_matrix_type Mpop (nb_indiv, 3); //Matrix of a solution

	//Array of matrix = list of the solutions
	Marray_type sol_list;

	int nb_indiv_pop = Mpop.size1();
    int nb_snp_pop = Mpop.size2();
    int indiv_select;


    std::cout << "nb indiv pop " << nb_indiv_pop << std::endl;
    std::cout << "nb snp pop " << nb_snp_pop << std::endl;

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
Marray_type cross_sol(Marray_type pop_sol){

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
	int_matrix_type MParent1 = pop_sol[parent1];
	std::cout << "parent 1 " << ":" << MParent1 << std::endl;
	int_matrix_type MParent2 = pop_sol[parent2];
	std::cout << "parent 2 " << ":" << MParent2 << std::endl;

	//second pair of parents
	int_matrix_type MParent3 = pop_sol[parent3];
	std::cout << "parent 3 " << ":" << MParent3 << std::endl;
	int_matrix_type MParent4 = pop_sol[parent4];
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

	int_matrix_type Mchild1;
	int_matrix_type Mchild2;
	int_matrix_type Mchild3;
	int_matrix_type Mchild4;

	for (int i = 0; i < int(MParent1.size1()); i++){
		for (int j = 0; j < cut1; j++){
			Mchild1(i, j) = MParent1(i, j);
		}
	}

	std::cout << "child 1  " << ":" << Mchild1 << std::endl;


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

    for (int i = 1; i < 20; i++)
    	std::cout << "Solution " << i << ":" << Apop[i] << std::endl;

    cross_sol(Apop);


}
