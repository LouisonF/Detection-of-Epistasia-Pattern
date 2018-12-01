/*
 * Miscellaneous.hpp
 *
 *  Created on: 24 nov. 2018
 *      Author: louison
 */

#ifndef MISCELLANEOUS_HPP_
#define MISCELLANEOUS_HPP_

#include <fstream>
#include <iostream>
#include <string>
#include <list>
#include <random>
#include <algorithm>

using namespace std;

class Miscellaneous {
public:
	Miscellaneous();
	virtual ~Miscellaneous();
	static void append_list(list<unsigned>, list<unsigned>);
	static void remove_list_from_list(list<unsigned>, list<unsigned>);
	static void append_to_file(string, string);
	static vector<vector<unsigned int>> combinator(vector<unsigned int> snps_sorted, unsigned int size);
	static void random_subset(vector<unsigned int> &in_subset, vector<unsigned int> &out_subset, unsigned int n_to_draw, mt19937 rand_seed);
	static vector<vector<unsigned int>> link_comb_to_snp(vector<unsigned int> snp_table, vector<vector<unsigned int>> all_index_combinations, unsigned int size);
	//vector<vector<unsigned int>> generate_all_combinations(vector<unsigned int> &snps_sorted,  int size);
	//void combinator(vector<vector<unsigned int>> output, vector<unsigned int> random_snp,vector<unsigned int> temp_combination, unsigned int i , int size);


};

#endif /* MISCELLANEOUS_HPP_ */
