/*
 * G2test.cpp
 *
 *  Created on: 23 nov. 2018
 *      Author: courtin
 */

#include "G2test.h"
#include <boost/math/distributions/chi_squared.hpp>

G2test::G2test(int_matrix_type cont_table, int_matrix_type theo_table) : cont_table(cont_table), theo_table(theo_table) {
	df = 0;
	pval = 0;
	g2 = 0;
}

G2test::~G2test() {
	// TODO Auto-generated destructor stub
}

void G2test::run_G2(){

	//    g2 = 0;      // initialized in constructor and this function is called in constructor
	for(unsigned i=0; i<cont_table.size1(); ++i)
	{
		for(unsigned j=0; j<cont_table.size2(); ++j)
		{
			if(cont_table(i,j) != 0 /*and cont_table(i,j) >= 5*/)
			{
				double div = (double) cont_table(i,j) / theo_table(i,j);
				g2 += cont_table(i,j) * log(div);
			}
		}
	}
	g2 *= 2;
	df = (cont_table.size1()-1)*(cont_table.size2()-1);
	boost::math::chi_squared_distribution<double> chi2_dist(df);
	pval = 1 - boost::math::cdf(chi2_dist, g2);
}

void G2test::display_g2(){
	cout << "g2 :" << g2 << endl;
	cout << "pval : " << pval << endl;
}

float G2test::get_g2(){
	return (g2);
}

float G2test::get_pval(){
	return (pval);
}
