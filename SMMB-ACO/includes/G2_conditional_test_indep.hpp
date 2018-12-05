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
	G2_conditional_test_indep(blas::matrix<int> const& genos,blas::matrix<int> const& phenos, vector<unsigned int> const& cond_genos);

    G2_conditional_test_indep(blas_column const& genos, blas_column const& phenos, blas_matrix const& cond_genos_vector);


    G2_conditional_test_indep(blas::matrix<int> genos,
    		blas::matrix<int> phenos,
                             vector<unsigned int> const& cond_genos_indexes, bool print_contingency=false);

    G2_conditional_test_indep(Contingency const& c, unsigned number_of_sub_contingencies);

    void run(bool verbose=false);
    void print_contingencies();

protected:
    std::vector<Contingency> _contingencies;
};

#endif // G2_CONDITIONAL_TEST_INDEP_HPP
