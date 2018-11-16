/*
 * Parametersfileparsing.cpp
 *
 *  Created on: 16 nov. 2018
 *      Author: Louison Fresnais, François Courtin
 *      Most of this code is from the SMMB-ACO implementation from Clément Niel
 */

#include "Parametersfileparsing.hpp"

#include <iostream>
#include <fstream>

using namespace std;
Parameters_file_parsing::Parameters_file_parsing(const string file_path)
{
	string file_path = file_path;

}

Parameters_file_parsing::~Parameters_file_parsing()
{
	// TODO Auto-generated destructor stub
}

//=================================================
// Parameters_file_parsing : Parsing
//=================================================


void Parsing()
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

void Parameters_file_parsing::import_line(string const& line)
{
    vector<string> token = this->split(line, ' ');
    string const& key = token[0];
    string & value = token[1];

    if(key == "header")
        header = atoi(value.c_str());

    else if(key == "separator")
    {
        if(value == "\t")
            separator = '\t';
        else
            separator = value.at(0);
    }
    else if(key == "alpha")

        alpha = atof(value.c_str());

    else if(key == "precision")
        precision = atof(value.c_str());

    else if(key == "subset_size_small")
        subset_size_small = atoi (value.c_str());

    else if(key == "n_trials_to_learn_mbs")
        n_trials_to_learn_mbs = atoi(value.c_str());

    else if(key == "n_trials_to_learn_1_mb")
        n_trials_to_learn_1_mb = atoi(value.c_str());

    else if(key == "gfile")
        genos_file = value;

    else if(key == "pfile")
        phenos_file = value;

    else if(key == "n_smmb_aco_runs")
        n_smmb_aco_runs = atoi(value.c_str());

    else if(key == "aco_n_ants")
        aco_n_ants = atoi(value.c_str());

    else if(key == "aco_set_size")
        aco_set_size = atoi(value.c_str());

    else if(key == "aco_n_iterations")
        aco_n_iterations = atoi(value.c_str());

    else if(key == "aco_tau_init")
        aco_tau_init = atof(value.c_str());

    else if(key == "aco_rho")
        aco_rho = atof(value.c_str());

    else if(key == "aco_lambda")
        aco_lambda = atof(value.c_str());

    else if(key == "aco_eta")
        aco_eta = atof(value.c_str());

    else if(key == "aco_alpha")
        aco_alpha = atof(value.c_str());

    else if(key == "aco_beta")
        aco_beta = atof(value.c_str());

    else {}

    n_mbs = number_ants * number_aco_iter;
}

//=================================================
// Parameters_file_parsing : split
//=================================================
//Method from Clément Niel's code
vector<string> Parameters_file_parsing::split(string const& s, char delim)
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
void Parameters_file_parsing::list_parameters() const
{
    cout << "########### PARAMETERS ###########\n" << "header => " << header << endl
    << "separator => " << separator << endl
    << "alpha => " << alpha << endl
    << "precision => " << precision << endl
    << "number_snp_per_ant => " << number_snp_per_ant << endl
    << "subset_size_small => " << subset_size_small << endl
    << "n_trials_to_learn_mbs => " << n_trials_to_learn_mbs << endl
    << "n_trials_to_learn_1_mb => " << n_trials_to_learn_1_mb << endl
    << "#################################" << endl;
}

//=================================================
// Parameters_file_parsing : update_subset_size_large
//=================================================
//Method from Clément Niel's code
void Parameters_file_parsing::update_subset_size_large(unsigned const& n_genos)
{
    if(number_snp_per_ant == 0)
    	number_snp_per_ant = sqrt(n_genos);
}

