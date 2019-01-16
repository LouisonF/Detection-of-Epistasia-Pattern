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
													 blas::matrix<int> ref_matrix,
                                                     bool do_print_contingency)
{
    cout << "METHOD G2_conditional_test_indep" << endl;
    unsigned n_obs = ref_matrix.size1();
    unsigned n_snp = ref_matrix.size2();
    unsigned n_cond_genos = cond_genos_indexes.size();
    cout << "there is : "<< n_cond_genos<< "conditional snps"<<endl;
    unsigned n_contingencies = pow(3,n_cond_genos);
    cout << "the number of contigency table is : "<<n_contingencies<<endl;
    if (n_contingencies == 0)
    {
    	n_contingencies=1;
    }

    blas_column genos_column = blas_column(genos,0);
    //blas::matrix<unsigned int> genos_current = genos_column.data();
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
        cout << "size of ref_genos_matrix" << ref_genos_matrix.size1()<<endl;
        cout << "size of ref_genos_matrix2" << ref_genos_matrix.size2()<<endl;
       // int val = ref_genos_matrix.size2();
        /*for(int h=0;h<=val;h++)
        {
        	cout << "#" << ref_genos_matrix(0,h)<<"#";
        }*/
		//TODO check if the data method is working

        for(unsigned i=0; i<n_obs; ++i)
        {
            // Put the current observation in the correct contingency table
            unsigned contingency_index = 0;
            unsigned int j=0;
            //for(int j=0;j<cond_genos_indexes.size(); ++j)
            for(auto it = cond_genos_indexes.begin(); it != cond_genos_indexes.end();it++,j++)
            {
            	//cout << "it equals to   " << *it << "   i equals to   " << i<<endl;
            	//cout << "j equals to    " << j<<endl;
                contingency_index += pow(3,j) * ref_matrix(i,*it);
            }
            //cout <<"size of con_genos_indexes : "<<cond_genos_indexes.size()<<endl;

            //cout << "CONTINGENCY INDEX OUI "<< contingency_index<<endl;
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

    run(do_print_contingency);
    cout << "METHOD G2_conditional_test_indep finished" << endl;
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

bool G2_conditional_test_indep::get_reliable()
{
	return reliable;
}
