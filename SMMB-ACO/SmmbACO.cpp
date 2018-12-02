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
    srand(time(NULL));
    rand_seed.seed(std::time(0));

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
		//We clear the score map(key > value) for the new iteration
		scores.clear();

		compute_distrib_prob();
		compute_cumlative_dristrib_proba();

		//Parrallel computation of each ant

		#pragma omp parallel for
        for(unsigned int ant=0; ant<_params.number_ants; ant++)
        {
            // Define a SNP set which is going to be selected by the sampling function
            vector<unsigned int> snp_table;
            snp_sampling(snp_table);

            // learn MB from these SNPS
            vector<unsigned int> mb;

            // add candidate mb to _mbs

        }


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

float Smmb_ACO::pheromone_for_snp(float tau_for_snp, float eta_for_snp)
{
	float effective_proba;
	tau_for_snp = pow(tau_for_snp, _params.aco_alpha)*pow(eta_for_snp, _params.aco_beta);

	effective_proba = tau_for_snp/sum_of_tau;

	return effective_proba;
}
//This method compute the distribution of probability
void Smmb_ACO::compute_distrib_prob()
{
	for(int i=0; i<tau.size(); i++)
	{
		pdf.at(i) = pheromone_for_snp(tau.at(i),eta.at(i));
	}
}
//This method ad to the cumulated distribution of probabilities table, evey proba that is greater than 0
void Smmb_ACO::compute_cumlative_dristrib_proba()
{
	float memory;
	cumulated_distrib_prob.clear();
	for(int i=0; i<pdf.size(); i++)
	{
		if (pdf.at(i) >0)
		{
			memory += pdf.at(i);
			cumulated_distrib_prob[memory].push_back(i);
		}
	}
}
unsigned int Smmb_ACO::select_snp_in_distrib_prob(float prob)
{
	//it_after_prob tell the iterator to start after a key equal of greater than prob
	//it's like for (i=2; i<vector.size(); i++) we don't start at 0 but at lower_bound set at 2
	auto it_after_prob = cumulated_distrib_prob.lower_bound(prob);
	// it_forward refers to the first key after the lower_bound.
	auto it_forward = it_after_prob++;
	//We are going to use the capacity of a map : first is going to
	//refers to the key (proba) and second is going to refers to the snp indice.
	if (it_after_prob != cumulated_distrib_prob.end())
	{
		auto it = it_forward;
		if((it_forward-> first) != prob)//if current key not equal to proba
		{
			it = it_after_prob;
			//it take the value of the lower_bound because the next key is not equal to prob
		}

		//Check that there is only one snp for this proba
		if((it->second.size() == 1))
		{
			//We push the snp into the selected snp vector
			return it->second.at(0);
		}else
		{
			//we have to randomly select one of the SNPs pointed by this key
			int rand_snp = rand() % it->second.size() + 1;
			return it->second.at(rand_snp);
		}
	}else
	{
		return 1;
	}

}

void Smmb_ACO::snp_sampling(vector<unsigned int> &snp_table)
{
	//SNP already in the snp_table ? true or false
	bool snp_in_sample = false;
	unsigned int i;

	//define the size of the snp_table
	unsigned int snp_table_size = _params.aco_set_size;

	for(i=0; i<snp_table_size; i++)
	{
		while(snp_in_sample == false)//while no snp added to snp_table, draw a new snp and test
		{
			float proba = ((float) rand() / (RAND_MAX));
			unsigned int drawn_snp = select_snp_in_distrib_prob(proba);

			//if nothing found, the find method return the last position of the vector
			if(find(snp_table.begin(), snp_table.end(), drawn_snp) == snp_table.end())
			{
				//if drawn_snp is not in snp_table, we add it
				snp_in_sample = true;
				snp_table.push_back(drawn_snp);

			}else
			{
				snp_in_sample = false;
			}
		}

	}

}

void Smmb_ACO::learn_mb(vector<unsigned int> &mb, vector<unsigned int> &snp_table)
{
	int counter;
	vector<unsigned int> memory_mb; //This a the vector of all the trials to learn a mb during this void
	while(counter < _params.max_trials_learn_mb )
	{
		memory_mb = mb;
		forward_phase(mb,snp_table);
		//backward_phase(mb,snp_table);
		counter++;
	}

	//call forward
	/*
	 * échantillonnage sur snp_table, selon taille sous-ensemble
	 * tri: tri pour pouvoir sort ensuite
	 *
	 */
}
void Smmb_ACO::forward_phase(vector<unsigned int> &mb, vector<unsigned int> &snp_table)
{
	vector<unsigned int> random_snps;
	//New sampling phase
	Miscellaneous::random_subset(snp_table,random_snps,_params.smallest_subset_size ,rand_seed);
	//Sorting is important for the combination method.
	sort(random_snps.begin(),random_snps.end());

	//Generating all combinations for the snp list
	unsigned int size = 3; //TODO maximum size of a combination, ask if we need to add it in parameters file
	vector<vector<unsigned int>> all_snps_combinations;
	Miscellaneous::combinator(random_snps,all_snps_combinations,size);
	Miscellaneous::link_comb_to_snp(random_snps,all_snps_combinations);


}
	//call backward
	/*
	 * pour tout X element de la MB
	 * 	pour toute combinaison S non vide, element de la MB
	 * 		test G2 entre X et T conditionnellement à S
	 * 		si p-valeur > alpha
	 * 			On enlève X de la MB, break;
	 * 		finsi
	 *
	 * 		Il faut faire des listes
	 */
void Smmb_ACO::backward_pahse(vector<unsigned int> &mb, vector<unsigned int> &snp_table)
{

}

//somme des tau DONE

// calcul des distributions de proba et des proba de chaque snp DONE

// TODO calcul parallelle et boucle pour faire travailler chaque fourmi.

//TODO select SNP, fonction de sampling
// TODO learn Markov blanket
// TODO ajout des mB à la secltion de mbs
//TODO mise a jour de tau
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
