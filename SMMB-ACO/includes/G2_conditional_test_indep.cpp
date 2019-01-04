#include "G2_conditional_test_indep.hpp"

#include <boost/numeric/ublas/io.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <iostream>

using namespace std;

//----------------------------------------------------
// G2_conditional_test_indep : Constructor 3
//----------------------------------------------------
G2_conditional_test_indep::G2_conditional_test_indep(blas::matrix<int> genos, blas::matrix<int> phenos,
                                                     vector<unsigned int> const& cond_genos_indexes,
                                                     bool do_print_contingency)
{
//    cout << "METHOD G2_conditional_test_indep" << endl;
    unsigned n_obs = genos.size1();
    unsigned n_cond_genos = cond_genos_indexes.size();
    unsigned n_contingencies = pow(3, n_cond_genos);
    blas_column genos_column = blas_column(genos,1);
    blas_column phenos_column = blas_column(phenos,1);
    _df = 2;
    if(n_cond_genos != 0)
        _df *= 3*n_cond_genos;

    _contingencies = vector<Contingency>(n_contingencies);

    // Fill contingency table (one or multiple)
    if(!cond_genos_indexes.empty())
    {
        blas::matrix<unsigned int> ref_genos_matrix;
		ref_genos_matrix = genos_column.data(); // get matrix from a column
        cout << "size of ref_genos_matrix" << ref_genos_matrix.size1()<<endl;
        cout << "size of ref_genos_matrix2" << ref_genos_matrix.size2()<<endl;
       // int val = ref_genos_matrix.size2();
        cout << ref_genos_matrix<<endl;
        /*for(int h=0;h<=val;h++)
        {
        	cout << "#" << ref_genos_matrix(0,h)<<"#";
        }*/
		//TODO check if the data method is working
        for(unsigned i=0; i<n_obs; ++i)
        {
            // Put the current observation in the correct contingency table
            unsigned contingency_index = 0;
            unsigned j=0;
            //cout <<"size of con_genos_indexes : "<<cond_genos_indexes.size()<<endl;
            for(vector<unsigned int>::const_iterator it=cond_genos_indexes.begin(); it!=cond_genos_indexes.end(); ++it, ++j)
            {
            	//cout << "it equals to   " << *it << "   i equals to   " << i<<endl;
            	//cout << "j equals to    " << j<<endl;
                contingency_index += pow(3, j) * ref_genos_matrix(i, 0);
            }

            Contingency& c = _contingencies[0];
            unsigned cr = phenos(i,0);
            unsigned cc = genos(i,0);
            c(cr, cc) += 1;
        }
    }
    else
    {
    	//TODO Test to compilate wihtout adding 0 to the contingency matrix
        for(unsigned i=0; i<n_obs; ++i)
        {
            Contingency& c = _contingencies[0];
            unsigned cr = phenos(i,0);
            unsigned cc = genos(i,0);
            c(cr, cc) += 1;
        }
    }

    run(do_print_contingency);
//    cout << "METHOD G2_conditional_test_indep finished" << endl;
}

//----------------------------------------------------
// G2_conditional_test_indep : Constructor 4
//----------------------------------------------------
G2_conditional_test_indep::G2_conditional_test_indep(Contingency const& c, unsigned number_of_sub_contingencies)
{
    unsigned nrows = c.size1();
    unsigned n_cols_by_contingency = c.size2() / number_of_sub_contingencies; //TODO check the number of sub_contingency table computation
    cout << "n_cols_by_contingency = " << n_cols_by_contingency << endl;

    _contingencies = vector<Contingency>(number_of_sub_contingencies);
    for(unsigned i=0; i<number_of_sub_contingencies; i++)
    {
        Contingency & current_c = _contingencies[i];
        current_c.resize(nrows, n_cols_by_contingency);
        for(unsigned j=0; j<nrows; j++)
        {
            for(unsigned k=i*n_cols_by_contingency, l=0; k<(i+1)*n_cols_by_contingency;k++,l++)
            {
                current_c(j,l) = c(j,k);
            }
        }
    }
    _df = 3;
    run(true);
}

//-----------------------------------------
// G2_conditional_test_indep : run
//-----------------------------------------
void G2_conditional_test_indep::run(bool verbose)
{
    _g2 = 0;
    for(unsigned i=0; i<_contingencies.size(); ++i)
    {
        G2_test_indep g2(_contingencies[i]);
        if(!g2.is_reliable())//TODO check if the is.reliableb method is working
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

//-----------------------------------------
// G2_conditional_test_indep : print_contingencies
//-----------------------------------------
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
