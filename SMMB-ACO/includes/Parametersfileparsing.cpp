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

Parameters_file_parsing::Parameters_file_parsing()
{
	// TODO Auto-generated destructor stub
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
    		sep = "\t"; // TODO: Problem here, need invalid conversion from char* to char
    	}
    	else
    	{
    		sep = value.at(0);
    	}
    }else if(key == "alpha")
    {
    	alpha = atof(value.c_str());
    }
    else if(key == "number_snp_per_ant")
    {
    	number_snp_per_ant = atoi(value.c_str());
    }
    else if(key == "precision")
    {
    	precision = atof(value.c_str());
    }
    else if(key == "smallest_subset_size")
    {
    	smallest_subset_size = atoi(value.c_str());
    }
    else if(key == "max_trials_smmb")
    {
    	max_trials_smmb = atoi(value.c_str());
    }
    else if(key == "max_trials_learn_mb")
    {
    	max_trials_learn_mb = atoi(value.c_str());
    }
    else if(key == "number_smmbaco_runs")
    {
    	number_smmbaco_runs = atoi(value.c_str());
    }
    else if(key == "number_ants")
    {
    	number_ants = atoi(value.c_str());
    }
    else if(key == "aco_set_size")
    {
    	aco_set_size = atoi(value.c_str());
    }
    else if(key == "number_aco_iter")
    {
    	number_aco_iter = atoi(value.c_str());
    }
    else if(key == "aco_tau_init")
    {
    	aco_tau_init = atof(value.c_str());
    }
    else if(key == "aco_rho")
    {
    	aco_rho = atof(value.c_str());
    }
    else if(key == "aco_lambda")
    {
    	aco_lambda = atof(value.c_str());
    }
    else if(key == "aco_eta")
    {
    	aco_eta = atof(value.c_str());
    }
    else if(key == "aco_alpha")
    {
    	aco_alpha = atof(value.c_str());
    }
    else if(key == "aco_beta")
    {
    	aco_beta = atof(value.c_str());
    }
    else if(key == "genos_file")
    {
    	genos_file = value.c_str();
    }
    else if(key == "phenos_file")
    {
    	phenos_file = value.c_str();
    }
    else
    {
    	cout << "parameter unknown"<<endl;
    }

    n_mbs = number_ants * number_aco_iter;
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
    << "alpha => " << alpha << endl
    << "precision => " << precision << endl
    << "number_snp_per_ant => " << number_snp_per_ant << endl
    << "smallest_subset_size => " << smallest_subset_size << endl
    << "max_trials_smmb => " << max_trials_smmb << endl
    << "max_trials_learn_mb => " << max_trials_learn_mb << endl
	<< "genos_file => " << genos_file << endl
	<< "phenos_file => "<< phenos_file <<endl
	<< "number_smmbaco_runs =>" << number_smmbaco_runs <<endl
	<< "number_ants =>" << number_ants <<endl
	<< "aco_set_size =>" << aco_set_size <<endl
	<< "number_aco_iter =>" << number_aco_iter <<endl
	<< "aco_tau_init => " << aco_tau_init <<endl
	<< "aco_rho => " << aco_rho <<endl
	<< "aco_lambda => " << aco_lambda <<endl
	<< "aco_eta => " << aco_eta <<endl
	<< "aco_alpha => " << aco_alpha <<endl
	<< "aco_beta => " << aco_beta <<endl
    << "#################################" << endl;
}
