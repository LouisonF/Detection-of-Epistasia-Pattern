/*
 * G2test.cpp
 *
 *  Created on: 23 nov. 2018
 *      Author: Louison Fresnais, François Courtin
 *      Project: SMMB-ACO and Genetic Algorithm for epistasis detection
 *      Under the supervision of Christine Sinoquet(Nantes University)
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

/*
 * *************************************************************
 */

//Test if the contingency table is reliable, return true if it is, flase if it is not
bool G2test::reliable_test(int_matrix_type & c)
{
	for(unsigned i=0; i<c.size1(); ++i)
	{
	   for(unsigned j=0; j<c.size2(); ++j)
	   {
           if(c(i,j) < 5)
			   return false;
	   }
	}
	return true;
}

/*
 * *************************************************************
 */

//Run the G2 test (based on Clement Niel code)
void G2test::run_G2(int &not_reliable, bool child){

	for(unsigned i=0; i<cont_table.size1(); ++i)
	{
		for(unsigned j=0; j<cont_table.size2(); ++j)
		{
			if(cont_table(i,j) != 0)
			{
				float div = (float) cont_table(i,j) / theo_table(i,j);
				g2 += cont_table(i,j) * log(div);
			}
		}
	}
	g2 *= 2;
	df = (cont_table.size1()-1)*(cont_table.size2()-1);


	//If the G2 different of infinity, compute the p-value
	if (g2 != std::numeric_limits<double>::infinity()){
		boost::math::chi_squared_distribution<double> chi2_dist(df);
		pval = 1 - boost::math::cdf(chi2_dist, g2);
	}else {
		//If the G2 is = to infinity, set the G2 to 0 to have a p-value of 1 and set it as a bad solution
		g2 = 0;
		boost::math::chi_squared_distribution<double> chi2_dist(df);
		pval = 1 - boost::math::cdf(chi2_dist, g2);
	}

	//Incrementation of non reliable test count if the test is not reliable
	if (!reliable_test(cont_table) and child){
		//cout << "Running G2..." << endl;
		not_reliable ++;
	}else{
		//cout << "Running G2..." << endl;
	}
}

/*
 * *************************************************************
 */

void G2test::display_g2(){
	cout << "g2 :" << g2 << endl;
	cout << "pval : " << pval << endl;
}

/*
 * *************************************************************
 */

double G2test::get_g2(){
	return (g2);
}

/*
 * *************************************************************
 */

double G2test::get_pval(){
	return (pval);
}
