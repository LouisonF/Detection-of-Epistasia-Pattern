/*
 * Parametersfileparsing.cpp
 *
 *  Created on: 16 nov. 2018
 *      Author: Louison Fresnais, François Courtin M2BB
 *      Project: SMMB-ACO and Genetic Algorithm for epistasis detection
 *      Under the supervision of Christine Sinoquet(Nantes University)
 *      Most of this code is from the SMMB-ACO implementation from Clément Niel
 *  Modified on: 19 nov. 2018
 */

#include "Parametersfileparsing.hpp"

#include <sstream>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <string>

using namespace std;

Parameters_file_parsing::Parameters_file_parsing(const string file_path): file_path(file_path)
{
	Parsing();
}


Parameters_file_parsing::~Parameters_file_parsing()
{
	// TODO Auto-generated destructor stub
}

//=================================================
// Parameters_file_parsing : Parsing
//=================================================


void Parameters_file_parsing::Parsing()
{
    ifstream content(file_path);
	if(content.is_open())
	{
		string line = "";
		while (getline(content,line))
		{
            if (line.length() != 0 && line[0] != '#')
            {
            	import_line(line);
            }
        }
    }
    else
    {
        std::cerr << "Error while opening your parameters file !\n";
    }
}
//=================================================
// Parameters_file_parsing : import_line
//=================================================

//Method from Clément Niel's code

void Parameters_file_parsing::import_line(string const line)
{
    vector<string> token = this->split(line, ' ');
    string const key = token[0];
    string  value = token[1];


    if(key == "header_nrows")
    {
    	header_nrows = atoi(value.c_str());
    }
    else if(key == "sep")
    {
    	if(value =="\t")
    	{
    		cerr << "tabulations are not supported as separator, use coma " << endl;
    	}
    	else
    	{
    		sep = value.at(0);
    	}
    }
    else if(key == "population_size")
    {
    	len_pop = atoi(value.c_str());
    	cout << "len pop" << endl;
    }
    else if(key == "pattern_size")
    {
    	len_pattern = atoi(value.c_str());
    	cout << "len pattern" << endl;
    }
    else if(key == "number_of_parents_selected")
    {
    	nb_parents = atoi(value.c_str());
    	cout << "nb parents" << endl;
    }
    else if(key == "number_of_iterations")
    {
    	nb_it = atoi(value.c_str());
    	cout << "nb it" << endl;
    }
    else if(key == "alpha")
    {
    	alpha = atof(value.c_str());
    	cout << "alpha" << endl;
    }
    else if(key == "mutation_probability")
    {
    	P_mutation = atof(value.c_str());
    	cout << "mutation" << endl;
    }
    else if(key == "bad_solution_selection_probability")
    {
    	P_selection = atof(value.c_str());
    	cout << "selection" << endl;
    }
    else
    {
    	cout << "parameter unknown"<<endl;
    }
}

//=================================================
// Parameters_file_parsing : split
//=================================================
//Method from Clément Niel's code
vector<string> Parameters_file_parsing::split(string const s, char delim)
{
    stringstream ss(s);
    string item;
    vector<string> tokens;
    while (getline(ss, item, delim))
        tokens.push_back(item);
    return tokens;
}

//=================================================
// Parameters_file_parsing : list_parameters
//=================================================

//Method from Clément Niel's code
void Parameters_file_parsing::list_parameters()
{
    cout << "########### PARAMETERS ###########\n" << "header_nrows => " << header_nrows << endl
    << "separator => " << sep << endl
    << "population_size => " << len_pop << endl
	<< "pattern_size => " << len_pattern << endl
    << "number_of_parents_selected => " << nb_parents << endl
    << "number_of_iterations => " << nb_it << endl
	<< "alpha => " << alpha << endl
    << "mutation_probability => " << P_mutation << endl
    << "bad_solution_selection_probability => " << P_selection << endl
    << "#################################" << endl;
}
