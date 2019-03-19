/*
 * G2test.h
 *
 *  Created on: 23 nov. 2018
 *      Author: Louison Fresnais, François Courtin
 *      Project: SMMB-ACO and Genetic Algorithm for epistasis detection
 *      Under the supervision of Christine Sinoquet(Nantes University)
 */

#ifndef G2TEST_H_
#define G2TEST_H_

#include "TypeDef.h"
using namespace std;

class G2test {
private:
	int df;
	double pval;
	double g2;
	int_matrix_type cont_table;
	int_matrix_type theo_table;
public:
	G2test(int_matrix_type,int_matrix_type);
	virtual ~G2test();
	bool reliable_test(int_matrix_type & c);
	void run_G2(int &not_reliable, bool child);
	void display_g2();
	double get_g2();
	double get_pval();
};

#endif /* G2TEST_H_ */
