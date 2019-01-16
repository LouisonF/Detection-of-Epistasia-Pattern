/*
 * Output.h
 *
 *  Created on: 3 déc. 2018
 *      Author: courtin
 */

#ifndef OUTPUT_H_
#define OUTPUT_H_

#include <iostream>
#include "TypeDef.h"

using namespace std;

class Output {
private:
	float_matrix_type Mpop_geno;
	vector<string> header;
	int len_pattern;
	int len_pop;
	string filename;
	vector<vector<float>> list_pattern;
	//vector<vector<float>> list_sol;
	//vector<vector<float>> best_pattern_list;

public:
	Output(float_matrix_type Mpop_geno, vector<string> header, int len_pattern, int len_pop, string filename);
	virtual ~Output();
	int cut(vector<vector<float>> & vect, int start, int end);
	void quickSort(vector<vector<float>> & vect, int p, int r);
	void set_list_pattern();
	//void set_list_sol();
	void set_best_sol();
	void write_best_sol();
};

#endif /* OUTPUT_H_ */
