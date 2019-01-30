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
		unsigned int subset_max_indice = in_subset.size();
		//cout << "subset_max_indice" << subset_max_indice <<endl;
		int random_indice = rand() % subset_max_indice + 0;
		//cout << "random indice is" << random_indice;
		if(find(out_subset.begin(), out_subset.end(), in_subset.at(random_indice)) == out_subset.end())
		{
			out_subset.push_back(in_subset.at(random_indice));
			value_in_vector = true;
		}else
		{
			value_in_vector = false;
		}
	}
	value_in_vector = false;

}

void Miscellaneous::combinator(vector<unsigned int> snps_sorted, vector<vector<unsigned int>> &all_index_combinations, unsigned int size)
{
    int n = snps_sorted.size();
    for(int s=1; s<=size; s++)
    {
    	cout << "value for n" << n <<endl;
        std::vector<bool> v(n);
        std::fill(v.begin(), v.begin() + s, true);

       do {
    	   vector<unsigned int> index_combination;
           for (int i = 0; i < n; ++i)
           {
               if (v[i])
               {
                   index_combination.push_back((i+1));
               }
           }
           all_index_combinations.push_back(index_combination);
       } while (std::prev_permutation(v.begin(), v.end()));
    }
    //TODO comment following block to remove combination verbose
       for (int x=0; x<all_index_combinations.size(); x++)
       {
    	   for(int y=0; y<all_index_combinations[x].size();y++)
    	   {
    		   cout << all_index_combinations[x][y];
    	   }
    	   cout << "\n";
       }

    //}

    //}
}
void Miscellaneous::link_comb_to_snp(vector<unsigned int> snps_sorted, vector<vector<unsigned int>> &all_index_combinations) //Link snp_sorted with their index
//This method is useless if all snps are in the subset but is usefull to match a subset of size < all snps with the input matrix
{
    for (int x=0; x<all_index_combinations.size(); x++)
    {
 	   for(int y=0; y<all_index_combinations[x].size();y++)
 	   {
 		  unsigned int tmp_indx =  all_index_combinations[x][y];
 		  all_index_combinations[x][y] = snps_sorted.at(tmp_indx-1);
 	   }
    }
}
void Miscellaneous::print_human_readable_combinations(vector<vector<unsigned int>> all_index_combinations)
{
    cout << "[" << endl;
    for(unsigned i=0; i<all_index_combinations.size(); i++)
    {
        vector<unsigned> const& currentV = all_index_combinations[i];
        cout << "\t[ ";
        for(unsigned j=0; j < currentV.size(); j++)
            cout << currentV[j] << " ";
        cout << "]" << endl;
    }
    cout << "]" << endl;
}

void Miscellaneous::append_vector_to_list(list<unsigned> & l, vector<unsigned> const& v)
{
    for(unsigned const& i: v)
        l.push_back(i);
}

bool Miscellaneous::compareFunc(pair<vector<unsigned>, vector<double>> const& a, pair<vector<unsigned>, vector<double>> const& b)
{
    if (a.second[0] < b.second[0])
    {
        return true;
    }
    else
    {
        return false;
    }
}
