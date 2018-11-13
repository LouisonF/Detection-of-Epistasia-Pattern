#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/array.hpp>

using namespace boost::numeric::ublas;

typedef boost::numeric::ublas::matrix<int> int_matrix_type;
typedef boost::array<int_matrix_type, 20> Marray_type;

Marray_type init_pop(int_matrix_type Mgeno, int nb_sol){



	int nb_indiv = Mgeno.size1();
	int nb_snp = Mgeno.size2();

	int_matrix_type Mpop (nb_indiv, 3);

	Marray_type sol_list;

	int nb_indiv_pop = Mpop.size1();
    int nb_snp_pop = Mpop.size2();
    int indiv_select;


    std::cout << "nb indiv pop " << nb_indiv_pop << std::endl;
    std::cout << "nb snp pop " << nb_snp_pop << std::endl;

    for (int i = 0; i < nb_sol; i++){
    	for (int j = 0; j < nb_indiv_pop; j++){
        	indiv_select = rand()%(nb_indiv);
        	Mpop (j, 0) = Mgeno(indiv_select, rand()%nb_snp);
        	Mpop (j, 1) = Mgeno(indiv_select, rand()%nb_snp);
        	Mpop (j, 2) = Mgeno(indiv_select, rand()%nb_snp);
    	}
    	sol_list[i] = {Mpop};
    }

    return sol_list;

}


Marray_type cross_sol(Marray_type sol_select){

}

int main () {

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

}
