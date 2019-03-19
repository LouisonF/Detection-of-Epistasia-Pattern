/*
 * Output.h
 *
 *  Created on: 3 déc. 2018
 *      Author: Louison Fresnais, François Courtin
 *      Project: SMMB-ACO and Genetic Algorithm for epistasis detection
 *      Under the supervision of Christine Sinoquet(Nantes University)
 */

#ifndef OUTPUT_H_
#define OUTPUT_H_

#include <iostream>
#include "TypeDef.h"

using namespace std;

class Output {
private:
	double_matrix_type Mpop_geno;
	vector<string> header;
	int len_pattern;
	int len_pop;
	string filename;
	vector<vector<double>> list_pattern;
	double alpha;
	//vector<vector<double>> list_sol;
	//vector<vector<double>> best_pattern_list;

public:
	Output(double_matrix_type Mpop_geno, vector<string> header, int len_pattern, int len_pop, string filename, double alpha);
	virtual ~Output();
	int cut(vector<vector<double>> & vect, int start, int end);
	void quickSort(vector<vector<double>> & vect, int p, int r);
	void set_list_pattern();
	//void set_list_sol();
	void set_best_sol();
	void write_best_sol();
};

#endif /* OUTPUT_H_ */
