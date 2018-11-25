/*
 * SmmbACO.cpp
 *
 *  Created on: 24 nov. 2018
 *      Author: louison
 */

#include "SmmbACO.hpp"

#include "Parametersfileparsing.hpp" // This include does not work, because it comes from a different project ... TODO


using namespace std;

Smmb_ACO::Smmb_ACO()
{
	//TODO Comment  passer un nom de fichier de paramètre sans le rentrer en dur ? Dans le script d'appel ?
	file_path = 'SMMB-ACO-parameters.txt';
	Parameters_file_parsing params(file_path);
	params.Parsing();
	params.list_parameters();
}

Smmb_ACO::~Smmb_ACO() {
	// TODO Auto-generated destructor stub
}

/*Tentative de définition de la marche à suivre ... */

void Smmb_ACO::run_aco()
{
	/*for(unsigned number_aco_it=0; number_aco_it<params.number_aco_iter; number_aco_it++)
	{
		cout << number_aco_it; //TODO: params instance is not declared in this scope. error in the constructor apparently.

	}*/
}

//BOUCLE FOR QUI FAIT TOURNER LA SMMB AVEC NMAX = LE NOMBRE D ITERATION

	//somme des tau

	// calcul des distributions de proba et des proba de chaque snp

	// calcul parallelle et boucle pour faire travailler chaque fourmis.

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
