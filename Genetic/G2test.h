/*
 * G2test.h
 *
 *  Created on: 23 nov. 2018
 *      Author: courtin
 */

#ifndef G2TEST_H_
#define G2TEST_H_

#include "TypeDef.h"
using namespace std;

class G2test {
private:
	int df;
	float pval;
	float g2;
	int_matrix_type cont_table;
	int_matrix_type theo_table;
public:
	G2test(int_matrix_type,int_matrix_type);
	virtual ~G2test();
	bool reliable_test(int_matrix_type & c);
	void run_G2(int &not_reliable, bool child);
	void display_g2();
	float get_g2();
	float get_pval();
};

#endif /* G2TEST_H_ */
