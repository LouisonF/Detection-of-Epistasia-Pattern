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

void Miscellaneous::combinator(vector<unsigned int> snps_sorted, vector<vector<unsigned int>> &all_index_combinations, unsigned int size)
{
    for(int s=1; s<=size; s++)
    {
    	int n = snps_sorted.size();
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
       /*for (int x=0; x<all_index_combinations.size(); x++)
       {
    	   for(int y=0; y<all_index_combinations[x].size();y++)
    	   {
    		   cout << all_index_combinations[x][y];
    	   }
    	   cout << "\n";
       }*/

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
/*
blas::matrix<float> Miscellaneous::init_contingency_table(int row,int col)
{

	blas::matrix<float> contingency_table(row,col);
	//Here we initialize the matrix with 0 in every cell
	for (int i = 0; i < int(contingency_table.size1()); i++){
		for (int j = 0; j < int(contingency_table.size2()); j++){
			contingency_table(i, j) = 0;
		}
	}
}
vector<blas::matrix<float>> Miscellaneous::contingency_table(vector<unsigned int> conditionnal_snp, blas::matrix<int> _genotypes, blas::matrix<int> _phenotypes, unsigned int number_of_obs_subset)
{
	vector<blas::matrix<float>> contingencies_vector;
	if(!conditionnal_snp.empty())
	{
		blas::matrix<int> temp
	}
}*/
/*DEPRECATED*/
	//Now we count the occurence for each possibility
/*	for (int i = 0; i < int(phenotypes.size1()); i++){
		for (int j = 0; j < int(genotypes.size2()); j++){
			if (phenotypes(i, 0) == 0){
				if(genotypes(i, j) == 0){
					contingency_table(0, 0) += 1;
				}
				if(genotypes(i, j) == 1){
					contingency_table(0, 1) += 1;
				}
				if(genotypes(i, j) == 2){
					contingency_table(0, 2) += 1;
				}
			}else{
				if(genotypes(i, j) == 0){
					contingency_table(1, 0) += 1;
				}
				if(genotypes(i, j) == 1){
					contingency_table(1, 1) += 1;
				}
				if(genotypes(i, j) == 2){
					contingency_table(1, 2) += 1;
				}
			}
		}
	}
	return(contingency_table);
}*/

//Sum of a matrix's row
/*unsigned int Miscellaneous::sum_row(blas::matrix<float> matrix, int row)
{
	float sum = 0;
	for (int i = 0; i < int(matrix.size2()); i++){
		sum += matrix(row, i);
	}
	return(sum);
}

//sum of a matrix's column
unsigned int Miscellaneous::sum_col(blas::matrix<float> matrix, int col)
{
	float sum = 0;
	for (int i = 0; i < int(matrix.size1()); i++){
		sum += matrix(i, col);
	}
	return(sum);
}

*/
//Calcul of a theorical contingency table from a contingency table
// i,j = (sum row i * sum col j)/total
/*blas::matrix<float> Miscellaneous::theorical_table(blas::matrix<float> table_cont)
{

	blas::matrix<float> table_theo(2, 3);
	unsigned int total = sum_row(table_cont, 0) + sum_row(table_cont, 1);
	for (int i = 0; i < int(table_theo.size1()); i++){
		for (int j = 0; j < int(table_theo.size2()); j++){
			table_theo(i, j) = (sum_row(table_cont, i) * sum_col(table_cont, j))/total;
		}
	}
	return(table_theo);
}*/
