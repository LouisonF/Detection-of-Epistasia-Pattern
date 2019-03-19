/*
 * Output.cpp
 *
 *  Created on: 3 déc. 2018
 *      Author: Louison Fresnais, François Courtin
 *      Project: SMMB-ACO and Genetic Algorithm for epistasis detection
 *      Under the supervision of Christine Sinoquet(Nantes University)
 */

#include "includes/Output.h"
#include <fstream>

Output::Output(double_matrix_type Mpop_geno, vector<string>header, int len_pattern, int len_pop, string filename, double alpha) : Mpop_geno(Mpop_geno), header(header), len_pattern(len_pattern), len_pop(len_pop), filename(filename), alpha(alpha) {
	// TODO Auto-generated constructor stub

}

Output::~Output() {
	// TODO Auto-generated destructor stub
}

/*
 * *************************************************************
 */

//Function for the Quicksort recursive function
int Output::cut(vector<vector<double>> & vect, int start, int end) {
	double pivot = vect[start][len_pattern+1];
	vector<double> vect_pivot(len_pattern+2);
	unsigned int from_left = start+1;
	unsigned int from_right = end;
	double tmp;

	for (int i = 0; i < len_pattern+2; i++){
		vect_pivot[i] = vect[start][i];
	}
	//cout << "Vector entering partition:";
	//print (a,start,end);
	//cout << endl;

	while (from_left != from_right) {
		if (vect[from_left][len_pattern+1]  <= pivot) from_left++;
		else {
			while (( from_left != from_right)  && (pivot < vect[from_right][len_pattern+1])) from_right--;
			for (int i = 0; i < len_pattern+2; i++){
				tmp =  vect[from_right][i];
				vect[from_right][i] = vect[from_left][i];
				vect[from_left][i] = tmp;

			}
		}
	}

	if (vect[from_left][len_pattern+1]>pivot) from_left--;
	for (int i = 0; i < len_pattern+2; i++){
		vect[start][i] = vect[from_left][i];
	}

	for (int i = 0; i < len_pattern+2; i++){
		vect[from_left][i] = vect_pivot[i];
	}
	//vect[from_left][len_pattern+1] = pivot;
	//cout << "Vector after partition:   ";
	//print (a,start,end);
	//cout << endl;

	return (from_left);
}

/*
 * *************************************************************
 */

//QuickSort recursive function
void Output::quickSort(vector<vector<double>> & vect, int p, int r) {
  if (p < r) {
    int q = cut(vect, p, r);
    quickSort(vect, p, q - 1);
    quickSort(vect, q + 1, r);
  }
}

/*
 * *************************************************************
 */

//Set a list of one copy of all the pattern contained in the population
void Output::set_list_pattern(){
	//cout << "Set list pattern for output..." << endl;
	vector<double> pattern;
	for (int i = 0; i < len_pop; i++){
		pattern.clear();
		for (int j = 0; j < len_pattern+2; j++){
			pattern.push_back(Mpop_geno(i,j));
		}
		//If the pattern isn't already in the list, add it to the list
		if ((count(list_pattern.begin(), list_pattern.end(), pattern) == 0) and (pattern[len_pattern+1] <= alpha)){
			list_pattern.push_back(pattern);
		}
	}
}

/*
 * *************************************************************
 */

//Call the quicksort function and print the solutions with their g2 score and p-values
void Output::set_best_sol(){
	/*vector<int> occurences_list;
	vector<double> best_pattern;
	int max_occ = 0;
	int cpt;
	for (int i = 0; i < int(list_pattern.size()); i++){
		cpt = count(list_sol.begin(), list_sol.end(), list_pattern[i]);
		occurences_list.push_back(cpt);
		if (cpt > max_occ){
			best_pattern_list.clear();
			max_occ = cpt;
			best_pattern = list_pattern[i];
			best_pattern_list.push_back(best_pattern);
		}
		else if (cpt == max_occ){
			best_pattern_list.push_back(list_pattern[i]);
		}

		cout << "pattern : ";
		for (int j = 0; j < len_pattern; j++){
			cout << list_pattern[i][j] << "/";
		}
		cout << " : " << occurences_list[i] << " occurence(s) " << "score : " << list_pattern[i][len_pattern] << " P_value : " << list_pattern[i][len_pattern+1] << endl;
	}
	cout << "Best pattern : " << endl;
	for (int i = 0;i < int(best_pattern_list.size()); i++){
		for (int j = 0; j < len_pattern; j++){
			cout << best_pattern_list[i][j] << "/";
		}
		cout << " with " << max_occ << " occurences, " << "score : " << best_pattern_list[i][len_pattern] << " P_value : " << best_pattern_list[i][len_pattern+1] << endl;
	}*/

	set_list_pattern();
	//cout << "Quick sorting..." << endl;
	quickSort(list_pattern, 0, list_pattern.size()-1);

	//cout << "List of solutions : " << endl;
	for (int i = 0; i < int(list_pattern.size()); i++){
		for (int j = 0; j < len_pattern; j++){
			//cout << list_pattern[i][j] << "/";
		}
		//cout << "\t score : " << list_pattern[i][len_pattern] << "\t p-value : " << list_pattern[i][len_pattern+1] << endl;
	}
}

/*
void Output::write_best_sol(){
	cout << "Write solutions..." << endl;

	ofstream file("/home/courtin/Documents/M2/ProjetC/detection-of-epistasia-pattern/Genetic_results/" + filename + ".txt", ios::out);

	if(file)
	{
		for (int i = 0;i < int(list_pattern.size()); i++){
			file << "{";
			for (int j = 0; j < len_pattern+2; j++){
				if (j == len_pattern -1 ){
					file << header[list_pattern[i][j]] << "}\t";
				}else if (j >= len_pattern){
					file << list_pattern[i][j] << "\t";
				}else{
					file << header[list_pattern[i][j]] << ",";
				}
			}
			file << endl;
		}


		file.close();
	}
	else
		cerr << "Cannot open file." << endl;
}
*/


/*
 * *************************************************************
 */

//Write the solutions, g2 score and p-value in the output file
void Output::write_best_sol(){

	//cout << "Write solutions..." << endl;
	ofstream file(filename + ".txt", ios::out);

	if(file)
	{
		for (int i = 0;i < int(list_pattern.size()); i++){
			file << "{";
			for (int j = 0; j < len_pattern+2; j++){
				if (j == len_pattern -1 ){
					file << list_pattern[i][j] << "}\t";
				}else if (j >= len_pattern){
					file << list_pattern[i][j] << "\t";
				}else{
					file << list_pattern[i][j] << ",";
				}
			}
			file << endl;
		}


		file.close();
	}
	else
		cerr << "Cannot open file." << endl;
}
