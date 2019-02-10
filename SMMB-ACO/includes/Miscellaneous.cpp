/*
 * Miscellaneous.cpp
 *
 *  Created on: 9 nov. 2018
 *      Author: Louison Fresnais, François Courtin
 *      Project: SMMB-ACO and Genetic Algorithm for epistasis detection
 *      Under the supervision of Christine Sinoquet(Nantes University)
 *  Modified on: 05 fev 2018
 */

#include "Miscellaneous.hpp"
using namespace std;

Miscellaneous::Miscellaneous() {
	// TODO Auto-generated constructor stub

}

Miscellaneous::~Miscellaneous() {
	// TODO Auto-generated destructor stub
}

/*
 *The random_subset method will take a vector of snps in input and output snps randomly selected in a vector of a given size.
 */

void Miscellaneous::random_subset(vector<unsigned int> &in_subset, vector<unsigned int> &out_subset, unsigned int n_to_draw, mt19937 rand_seed)
{
	//We generate a discrete distribution of in_subset size
	//Then we randomly, according to the seed, select a number in this distribution
	bool value_in_vector = false;
	//Generate an uniform distribution starting from 0 to the size of the input subset minus one
	uniform_int_distribution<unsigned int> dis(0, in_subset.size()-1);
	//While the output vector haven't reach the required size, random pick in input vector.
	while(out_subset.size() < n_to_draw)
	{
		unsigned int subset_max_indice = in_subset.size();
		//Random pick an int in the range(begin to subset_max_indice) of the input vector. this int is a snp indice.
		int random_indice = rand() % subset_max_indice + 0;
		//If the picked snp is not already in the vector, we push it in the output vector and indicate that the value is in the vector.
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

/*
 *The random_subset method will find all the combinations available for a epistasis pattern of a given size.
 */

void Miscellaneous::combinator(vector<unsigned int> snps_sorted, vector<vector<unsigned int>> &all_index_combinations, unsigned int size)
{
	//We get the number of snps in the previously radnomly picked snp vector
    int n = snps_sorted.size();
    //for each size of pattern 1 and 2 for exemple, seek all possible combinations.
    for(unsigned int s=1; s<=size; s++)
    {
    	//We create a vector of booleans of size n, the number of snps in the vector. We fill it with true.
        std::vector<bool> v(n);
        std::fill(v.begin(), v.begin() + s, true);
        //we iterate over the boolean vector and push the combination to a temporary index_combination vector.
        //When all combinations have been found for a size of pattern we push the index_combination vector to the all_index_combinations
        //vector of vectors.
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
       //while all permutations haven't been found according to the boolean vector size, iterate.
    }
    //TODO comment following block to remove combination verbose
       /*for (unsigned int x=0; x<all_index_combinations.size(); x++)
       {
    	   for(unsigned int y=0; y<all_index_combinations[x].size();y++)
    	   {
    		   cout << all_index_combinations[x][y];
    	   }
    	   cout << "\n";
       }*/
}

/*
 *The link_comb_to_snp method we link the combinations from combinator to the relevant indices in the snps_sorted vector.
 */

void Miscellaneous::link_comb_to_snp(vector<unsigned int> snps_sorted, vector<vector<unsigned int>> &all_index_combinations)
{
	//Get the values of the vector of vectors.
    for (int x=0; x<all_index_combinations.size(); x++)
    {
 	   for(int y=0; y<all_index_combinations[x].size();y++)
 	   {
 		   //replace the value in all_index_combinations by the snp indice value.
 		  unsigned int tmp_indx =  all_index_combinations[x][y];
 		  all_index_combinations[x][y] = snps_sorted.at(tmp_indx-1);
 	   }
    }
}

/*
 *The print_human_readable_combinations method can be called in order to see the combinations that are made by the algorithm.
 *Not used in normal usage of the algorithm.
 */

void Miscellaneous::print_human_readable_combinations(vector<vector<unsigned int>> all_index_combinations)
{
    cout << "[" << endl;
    //iterate over the all_index_combination vector's size.
    for(unsigned i=0; i<all_index_combinations.size(); i++)
    {
        vector<unsigned> const& currentV = all_index_combinations[i];
        cout << "\t[ ";
        //print the content of the current combination
        for(unsigned j=0; j < currentV.size(); j++)
            cout << currentV[j] << " ";
        cout << "]" << endl;
    }
    cout << "]" << endl;
}

/*
 *The compareFunc function is used to sort the map. THis function takes two pair of variables that refers to a map key/value.
 *These pairs of variable are then compared by their p-value (first value of the values vector in the results map).
 */

bool Miscellaneous::compareFunc(pair<vector<unsigned>, vector<double>> const& a, pair<vector<unsigned>, vector<double>> const& b)
{
	//if p-value of this mb is smaller than the p-value of the other mb, return true. else return false
	//(second mb have a smaller p-value than first)
    if (a.second[0] < b.second[0])
    {
        return true;
    }
    else
    {
        return false;
    }
}
