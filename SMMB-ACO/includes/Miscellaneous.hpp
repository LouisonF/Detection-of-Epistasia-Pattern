/*
 * Miscellaneous.hpp
 *
 *  Created on: 9 nov. 2018
 *      Author: Louison Fresnais, François Courtin
 *      Project: SMMB-ACO and Genetic Algorithm for epistasis detection
 *      Under the supervision of Christine Sinoquet(Nantes University)
 *  Modified on: 05 fev 2018
 */

#ifndef MISCELLANEOUS_HPP_
#define MISCELLANEOUS_HPP_

#include <fstream>
#include <iostream>
#include <string>
#include <list>
#include <random>
#include <algorithm>
#include "datainput.hpp"

using namespace std;

class Miscellaneous {
public:
	Miscellaneous();
	virtual ~Miscellaneous();
	static void combinator(vector<unsigned int> snps_sorted, vector<vector<unsigned int>> &all_index_combinations, unsigned int size);
	static void random_subset(vector<unsigned int> &in_subset, vector<unsigned int> &out_subset, unsigned int n_to_draw, mt19937 rand_seed);
	static void link_comb_to_snp(vector<unsigned int> snps_sorted, vector<vector<unsigned int>> &all_index_combinations); //Link snp_sorted with their index
	static void print_human_readable_combinations(vector<vector<unsigned int>> all_index_combinations);
	static bool compareFunc(pair<vector<unsigned>, vector<double>> const& a, pair<vector<unsigned>, vector<double>> const& b);
};

#endif /* MISCELLANEOUS_HPP_ */
