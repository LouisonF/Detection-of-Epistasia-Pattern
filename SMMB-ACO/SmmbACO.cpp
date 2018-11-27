/*
 * SmmbACO.cpp
 *
 *  Created on: 24 nov. 2018
 *      Author: louison
 */

#include "SmmbACO.hpp"
using namespace std;

//Constructor
Smmb_ACO::Smmb_ACO(blas::matrix<int> & genos, blas::matrix<int> & phenos, Parameters_file_parsing params): _genotypes(genos), _phenotypes(phenos), _params(params)
{

	//Variable assignation
	int number_of_snps;
	//Output file opening
    string file_basename = basename((char*)_params.genos_file.c_str());
    string result_filename = "outputs/RESULT_" + file_basename;
    output_file.open(result_filename.c_str(), ios::trunc);
    cout << "file_basename : "<< file_basename << "result_filename : " << result_filename << endl;

    if(!output_file)
    {
        std::cerr << "Error while opening output.txt (by writing access) !\n";
        exit(-1);
    }

	//Random seed
	srand (time(NULL));
	rand_seed = rand(); // TODO: Check if this is a good method to generate a seed.

    // Here, we count the number of independant test made during a SMMB execution.
    number_of_indep_test = 0;
    //sum_of_tau refers to the sum of pheromones levels
    sum_of_tau= 0.0;

    //number_of_snps refers to the number of snps available in the dataset
    number_of_snps = _genotypes.size2();

    tau.assign(number_of_snps,_params.aco_tau_init);
    eta.assign(number_of_snps,_params.aco_eta);
    pdf.assign(number_of_snps, 0.0);
    number_executions = params.number_smmbaco_runs;

}

//Destructor
Smmb_ACO::~Smmb_ACO()
{
	// TODO Auto-generated destructor stub
}

/*Tentative de définition de la marche à suivre ... */

void Smmb_ACO::run_ACO()
{
	cout << "SMMB is currently running" << endl;
	for(unsigned number_aco_it=0; number_aco_it< _params.number_aco_iter; number_aco_it++)
	 {
		cout << "ACO iteration number: " << number_aco_it << endl;
	 //Computation of the probability distribution
		//Sum of pheromones
		sum_tau();

	 }
}

//BOUCLE FOR QUI FAIT TOURNER LA SMMB AVEC NMAX = LE NOMBRE D ITERATION
void Smmb_ACO::sum_tau()
{
    // initialisation of the variable sum_of_tau at 0 (this void is called at each ACO iteration)
	sum_of_tau= 0.0;
	//Tried to use list type for tau and eta variable but list don't use direct access, vector does.
    for(auto i=0; i<tau.size(); i++)
    {
    	sum_of_tau += pow(tau.at(i), _params.aco_alpha) * pow(eta.at(i), _params.aco_beta);
    }

}
//somme des tau

// calcul des distributions de proba et des proba de chaque snp

// calcul parallelle et boucle pour faire travailler chaque fourmi.

//select SNP, fonction de sampling
//learn Markov blanket
//ajout des mv à la secltion de mbs
//mise a jour de tau
//Fin de la fonction qui fait tourner la smmb

//il faut une fonction de calcul du taux de phéromone.
//Une fonction de stockage en mémoire des proba et du taux de phréromone?
//Une fonction de sampling des snps: ceux qui n'ont pas encore été sample
//sont ajoutés à la selection de snps
//Une fonction qui permet de tirer les SNPs selon leur proba et au sort si deux snps ont la meme proba
//Une fonction qui permet l'update de tau (pas obligatoire mais allège le code)
//Une fonction qui permet la somme des taux de pheromone pour tous les sNPs
//Une fonction qui va concretement réaliser la markov blanket en appellant la phase forward et la phase backward
//Phase forward
//Phase backward
// Une fonction qui permet la raffinnement des résultats
// via la transformation de MB en une seule
// via la selection de SNPs uniques
// Réalisation des tests statistiques (G2)
