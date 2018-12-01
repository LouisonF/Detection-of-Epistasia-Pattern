/*
 * Miscellaneous.cpp
 *
 *  Created on: 24 nov. 2018
 *      Author: louison
 */

#include "Miscellaneous.hpp"
using namespace std;

Miscellaneous::Miscellaneous() {
	// TODO Auto-generated constructor stub

}

Miscellaneous::~Miscellaneous() {
	// TODO Auto-generated destructor stub
}

void Miscellaneous::append_list(list<unsigned> list1, list<unsigned> list2)
{
	list<unsigned> buffer = list2;
	list1.insert(list1.end(),buffer.begin(), buffer.end());
}

void Miscellaneous::remove_list_from_list(list<unsigned> list1, list<unsigned> list2)
{
	for(auto it=list1.begin(); it!=list1.end(); ++ it)
	{
		unsigned temp = *it;
		list2.remove(temp);
	}
}

void Miscellaneous::append_to_file(string filename, string data_to_append)
{
	ofstream content;
	content.open(filename, ofstream::out | ofstream::app);
	//option out allow to open the file in writing mode.
	//option app allow to append data at the end of the file.

	if(content.is_open())
	{
		content << data_to_append;
	}
	else
	{
		cerr << "error in file parsing \n";
		exit(-1);
	}
	content.close();
}

void Miscellaneous::random_subset(vector<unsigned int> &in_subset, vector<unsigned int> &out_subset, unsigned int n_to_draw, mt19937 rand_seed)
{
	//We generate a discrete distribution of in_subset size
	//Then we randomly, according to the seed, select a number in this distribution
	bool value_in_vector = false;
	uniform_int_distribution<unsigned int> dis(0, in_subset.size()-1);

	while(out_subset.size() < n_to_draw)
	{

		int random_indice = dis(rand_seed);
		if(find(out_subset.begin(), out_subset.end(), random_indice) == out_subset.end())
		{
			out_subset.push_back(in_subset.at(random_indice));
			value_in_vector = true;
		}else
		{
			value_in_vector = false;
		}
	}

}

vector<vector<unsigned int>> Miscellaneous::combinator(vector<unsigned int> snps_sorted, unsigned int size)
{
    for(int s=1; s<=size; s++)
    {
    	//cout << "size = " << s << endl;
    	vector<unsigned int> index_combination;
    	vector<vector<unsigned int>> all_index_combinations;
    	vector<bool> temp_vec(snps_sorted.size());
    	fill(temp_vec.end() - s, temp_vec.end(), true);

    	do {
    		for (int i = 0; i < snps_sorted.size(); ++i)
    		{
    			if (temp_vec[i])
    			{
    				cout << "temp_vec = " << temp_vec.at(i) << endl;
    				//cout << (i + 1) << " ";
    				index_combination.push_back(i+1);
    				cout << "index a i = " << index_combination.at(0) <<endl;
    			}
    		}
    		cout << "\n";
    	} while (next_permutation(temp_vec.begin(), temp_vec.end()));

    	all_index_combinations.push_back(index_combination);
    }
    return all_index_combinations;
}

vector<vector<unsigned int>> Miscellaneous::link_comb_to_snp(vector<unsigned int> snp_table, vector<vector<unsigned int>> all_index_combinations, unsigned int size)
{
	vector<vector<unsigned int>> all_snp_combinations;
    for(int s=1; s<=size; s++)
    {
    	//TODO : Faire le lien entre l'index d'un snp et sa valeur dans snp_table
    	//TODO : Comment grouper les combinaisons ...
    }
    return all_snp_combinations;
}
/*DEPRECATED
vector<vector<unsigned int>>  Miscellaneous::generate_all_combinations(vector<unsigned int> &snps_sorted, int size)
{
	cout << "ALL COMBINATIONS" <<endl;
	vector<vector<unsigned int>> all_combinations;
	vector<unsigned int> temp_combinations;
	unsigned int i = 0;
	combinator(all_combinations,snps_sorted,temp_combinations,i,size);
	return all_combinations;
}*/
/*DEPRECATED
void Miscellaneous::combinator(vector<vector<unsigned int>> output, vector<unsigned int> random_snp,vector<unsigned int> temp_combination, unsigned int i ,int size)
{
	cout << "COMBINATOOORRR" << endl;
	//if the combination is at the required size, we push it in the output vector.
	//size that -1 at each call of the recursive method
	if (size == 0)
	{
		output.push_back(temp_combination);
		cout << "condition de fin trouvée" << endl;
		return;
	}

	for (auto i=0; i != random_snp.size(); i++)
	{

		unsigned int tmp = random_snp.at(i);
		cout <<"temp = " << tmp << endl;
		cout << "i = " << i << endl;
		temp_combination.push_back(tmp);
		cout << "size = " << size << endl;


		combinator(output,random_snp,temp_combination,i+1,size-1);

		// Popping out last inserted element
		// from the vector
		temp_combination.pop_back();
	}
}*/

