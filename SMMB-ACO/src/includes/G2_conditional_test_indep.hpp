/*
 * G2_conditional_test_indep.hpp
 *
 *  Created on: 9 nov. 2018
 *      Author: Louison Fresnais, François Courtin
 *      Project: SMMB-ACO and Genetic Algorithm for epistasis detection
 *      Under the supervision of Christine Sinoquet(Nantes University)
 *  Modified on: 05 fev 2018
 */

#ifndef G2_CONDITIONAL_TEST_INDEP_HPP
#define G2_CONDITIONAL_TEST_INDEP_HPP

#include "G2_test_indep.hpp"
#include "Contingency.hpp"

#include <vector>
#include <list>

typedef blas::matrix<int> blas_matrix;
typedef blas::matrix_column<blas::matrix<int>> blas_column;
typedef blas::matrix_row<blas::matrix<int>> blas_row;

using namespace std;

class G2_conditional_test_indep : public G2_test_indep
{
public:

    G2_conditional_test_indep(blas::matrix<int> genos,
    		blas::matrix<int> phenos,
            vector<unsigned int> const& cond_genos_indexes,
			blas::matrix<int> ref_matrix,
			bool print_contingency=false);

    void run(bool verbose=false);
    void print_contingencies();
    bool get_reliable();

protected:
    std::vector<Contingency> _contingencies;
    bool reliable;
};

#endif // G2_CONDITIONAL_TEST_INDEP_HPP
