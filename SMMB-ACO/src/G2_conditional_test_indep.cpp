/*
 * G2_conditional_test_indep.cpp
 *
 *  Created on: 9 nov. 2018
 *      Author: Louison Fresnais, François Courtin
 *      Project: SMMB-ACO and Genetic Algorithm for epistasis detection
 *      Under the supervision of Christine Sinoquet(Nantes University)
 *  Modified on: 05 fev 2018
 */

//Reworked module from Clement Niel's code
#include "includes/G2_conditional_test_indep.hpp"
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <iostream>

using namespace std;

//----------------------------------------------------
// G2_conditional_test_indep : Constructor 1
//----------------------------------------------------
//This constructor call every required method to compute contingency tables and then run a g2_cond_test on contingency tables.
G2_conditional_test_indep::G2_conditional_test_indep(blas::matrix<int> genos, blas::matrix<int> phenos,
                                                     vector<unsigned int> const& cond_genos_indexes,
													 blas::matrix<int> ref_matrix,
                                                     bool do_print_contingency)
{

    unsigned n_obs = ref_matrix.size1();
    unsigned n_snp = ref_matrix.size2();
    unsigned n_cond_genos = cond_genos_indexes.size();
    unsigned n_contingencies = pow(3,n_cond_genos);
    if (n_contingencies == 0)
    {
    	n_contingencies=1;
    }
    blas_column genos_column = blas_column(genos,0);
    blas_column phenos_column = blas_column(phenos,0);
    _df = 2;
    if(n_cond_genos != 0)
        _df *= 3*n_cond_genos;

    _contingencies = vector<Contingency>(n_contingencies);

    // Fill contingency table (one or multiple)
    if(!cond_genos_indexes.empty())
    {
        blas::matrix<unsigned int> ref_genos_matrix;
		ref_genos_matrix = genos; // get matrix from a column
        for(unsigned i=0; i<n_obs; ++i)
        {
            // Put the current observation in the correct contingency table
            unsigned contingency_index = 0;
            unsigned int j=0;
            for(auto it = cond_genos_indexes.begin(); it != cond_genos_indexes.end();it++,j++)
            {
                contingency_index += pow(3,j) * ref_matrix(i,*it);
            }
            //Add +1 in to the relevant combination.
            Contingency& c = _contingencies[contingency_index];
            unsigned cr = phenos_column(i);
            unsigned cc = genos_column(i);
            c(cr, cc) += 1;
        }
    }
    else
    {
        for(unsigned i=0; i<n_obs; ++i)
        {
            Contingency& c = _contingencies[0];
            unsigned cr = phenos_column(i);
            unsigned cc = genos_column(i);
            c(cr, cc) += 1;
        }
    }
//run g2 on contingency tables
    run(do_print_contingency);
}
/*
 * *************************************************************
 */

//-----------------------------------------
// G2_conditional_test_indep : run
//-----------------------------------------
//run a g2_conditional_test_indep on a vector of contingency tables
void G2_conditional_test_indep::run(bool verbose)
{
    _g2 = 0;
    for(unsigned i=0; i<_contingencies.size(); ++i)
    {
        G2_test_indep g2(_contingencies[i]);
        if(!g2.is_reliable())
        {
            _g2 = 0;
            _pval = 1;
            _reliable = false;
            break;
        }
        _reliable = true;
        _g2 += g2.g2();
    }
    if(verbose)
    {
        cout << "g2_stat\t" << _g2 << "\ndf\t" << _df << endl;
        print_contingencies();
    }
    boost::math::chi_squared_distribution<double> chi2_dist(_df);
    _pval = 1 - boost::math::cdf(chi2_dist, _g2);
    if(_pval == 0)
        _pval = 2.0e-16;
}
/*
 * *************************************************************
 */

//-----------------------------------------
// G2_conditional_test_indep : print_contingencies
//-----------------------------------------
//if requested by the user in the g2_cond call, print contingencies tables.
void G2_conditional_test_indep::print_contingencies()
{
    cout << "Contingencies:" << endl;
    for(unsigned i=0; i<_contingencies.size(); i++)
    {
        Contingency& c = _contingencies[i];
        cout << c << endl;
    }
    cout << "end" << endl;
}
/*
 * *************************************************************
 */
 
bool G2_conditional_test_indep::get_reliable()
{
	return reliable;
}
