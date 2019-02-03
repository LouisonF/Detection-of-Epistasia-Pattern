/*
 * Output.cpp
 *
 *  Created on: 3 déc. 2018
 *      Author: courtin
 */

#include "Output.h"
#include <fstream>

Output::Output(float_matrix_type Mpop_geno, vector<string>header, int len_pattern, int len_pop, string filename) : Mpop_geno(Mpop_geno), header(header), len_pattern(len_pattern), len_pop(len_pop), filename(filename) {
	// TODO Auto-generated constructor stub

}

Output::~Output() {
	// TODO Auto-generated destructor stub
}

int Output::cut(vector<vector<float>> & vect, int start, int end) {
	float pivot = vect[start][len_pattern+1];
	unsigned int from_left = start+1;
	unsigned int from_right = end;
	float tmp;

	//cout << "Vector entering partition:";
	//print (a,start,end);
	//cout << endl;

	while (from_left != from_right) {
		if (vect[from_left][len_pattern+1]  <= pivot) from_left++;
		else {
			while (( from_left != from_right)  && (pivot < vect[from_right][len_pattern+1])) from_right--;
			//cout << "swaping " << a[from_left] << " with "<< a[from_right] << endl;
			tmp =  vect[from_right][len_pattern+1];
			vect[from_right][len_pattern+1] = vect[from_left][len_pattern+1];
			vect[from_left][len_pattern+1] = tmp;
		}
	}

	if (vect[from_left][len_pattern+1]>pivot) from_left--;
	vect[start][len_pattern+1] = vect[from_left][len_pattern+1];
	vect[from_left][len_pattern+1] = pivot;

	//cout << "Vector after partition:   ";
	//print (a,start,end);
	//cout << endl;

	return (from_left);
}


void Output::quickSort(vector<vector<float>> & vect, int p, int r) {
  if (p < r) {
    int q = cut(vect, p, r);
    quickSort(vect, p, q - 1);
    quickSort(vect, q + 1, r);
  }
}


void Output::set_list_pattern(){
	vector<float> pattern;
	for (int i = 0; i < len_pop; i++){
		pattern.clear();
		for (int j = 0; j < len_pattern+2; j++){
			pattern.push_back(Mpop_geno(i,j));
		}
		if (count(list_pattern.begin(), list_pattern.end(), pattern) == 0){
			list_pattern.push_back(pattern);
		}
	}
}

/*
void Output::set_list_sol(){

	vector<float> sol;
	for (int i = 0; i < len_pop; i++){
		sol.clear();
		for (int j = 0; j < len_pattern+2; j++){
			sol.push_back(Mpop_geno(i,j));
		}
		list_sol.push_back(sol);
	}
}
*/
void Output::set_best_sol(){
	/*vector<int> occurences_list;
	vector<float> best_pattern;
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
	quickSort(list_pattern, 0, list_pattern.size()-1);

	cout << "List of solutions : " << endl;
	for (int i = 0; i < list_pattern.size(); i++){
		for (int j = 0; j < len_pattern; j++){
			cout << list_pattern[i][j] << "/";
		}
		cout << "\t score : " << list_pattern[i][len_pattern] << "\t p-value : " << list_pattern[i][len_pattern+1] << endl;
	}
}

/*void Output::write_best_sol(){
	ofstream file("/home/courtin/M2/ProjetC/detection-of-epistasia-pattern/Genetic_results/" + filename + ".txt", ios::out);

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
}*/

void Output::write_best_sol(){
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
