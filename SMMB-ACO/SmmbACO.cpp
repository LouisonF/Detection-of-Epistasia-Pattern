/*
 * SmmbACO.cpp
 *
 *  Created on: 24 nov. 2018
 *      Author: louison
 */

#include "SmmbACO.hpp"
using namespace std;

//Constructor
Smmb_ACO::Smmb_ACO(blas::matrix<int> & genos, blas_column & phenos, Parameters_file_parsing params): _genotypes(genos), _phenotypes(phenos), _params(params)
{

	//Output file opening
    string file_basename = basename((char*)_params.genos_file.c_str());
    string result_filename = "outputs/RESULT_" + file_basename;
    //output_file.open(result_filename.c_str(), ios::trunc);
    cout << "file_basename : "<< file_basename << "result_filename : " << result_filename << endl;

    /*if(!output_file)
    {
        std::cerr << "Error while opening output.txt (by writing access) !\n";
        exit(-1);
    }*/

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

		//#pragma omp parallel for
        for(unsigned int ant=0; ant<_params.number_ants; ant++)
        {
            // Define a SNP set which is going to be selected by the sampling function
            vector<unsigned int> snp_table;
            snp_sampling(snp_table);

            // learn MB from these SNPS
            vector<unsigned int> mb;
            learn_mb(mb,snp_table);

            // add candidate mb to _mbs
            if(!mb.empty())
            {
            	mbs.push_back(mb);
            }
        }
        update_tau(); //TODO


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
void Smmb_ACO::evaporation_rate_update(unsigned int snp_index, float g2_score)
{
	tau[snp_index] = (1-_params.aco_rho) * tau[snp_index] + g2_score*_params.aco_lambda;
}

void Smmb_ACO::update_tau()
{
	for(auto it = scores.begin(); it != scores.end(); it++)
	{
		unsigned int snp_index = it->first;
		vector<float> scores_vec = it->second;
		for(unsigned int i = 0; i < scores_vec.size(); i++)
		{
			evaporation_rate_update(snp_index,scores_vec[i]);
		}
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
	unsigned int counter = 0;
	vector<unsigned int> memory_mb; //This a the vector of all the trials to learn a mb during this void
	while(counter < _params.max_trials_learn_mb )
	{
        cout << "hello, i am learning a markov blanket" <<endl; //TODO DEBUG
		memory_mb = mb;
		forward_phase(mb,snp_table);
		backward_phase(mb,snp_table);
		counter++;
	}
	backward_phase(mb,snp_table);
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
	float best_subset_pvalue;
	unsigned int best_snp_index;
	//New sampling phase
	Miscellaneous::random_subset(snp_table,random_snps,_params.smallest_subset_size ,rand_seed);
	//Sorting is important for the combination method.
	sort(random_snps.begin(),random_snps.end());

	//Generating all combinations for the snp list
	unsigned int size = 3; //TODO maximum size of a combination, ask if we need to add it in parameters file
	vector<vector<unsigned int>> all_snps_combinations;
	Miscellaneous::combinator(random_snps,all_snps_combinations,size);
	Miscellaneous::link_comb_to_snp(random_snps,all_snps_combinations);

	for (unsigned int i=0; i<all_snps_combinations.size(); i++)
	{
		vector<unsigned int> current_combination = all_snps_combinations[i];
		for(auto it =current_combination.begin(); it != current_combination.end(); it++)
		{
			vector<unsigned int> current_combination_temp = current_combination;
			unsigned int current_SNP = *it;
			cout << "current SNP is :" <<*it<<endl; //TODO DEBUG
			//We merge the temporary markov blanket
			vector<unsigned int> mb_temp = mb;
			//erase current_SNP from temporary combination vector
			current_combination_temp.erase(it);
			for(auto i=0; i<current_combination_temp.size(); i++)
			{
				mb_temp.push_back(current_combination_temp[i]);
			}

			//On fait le g2 conditionnel
			//On lui passe le genotype, le snp à évaluer, les phénotypes et la markov blaneket temporaire.
			// Conditional independance test
			//Getting the column
			blas::matrix_column<blas::matrix<int> > boostgenotype_column(blas_column(_genotypes,current_SNP));
			/*for(unsigned int i=0; i<_genotypes.size2();i++)
			{
				if(i == current_SNP)
				{
					for(unsigned int j=0; j<_genotypes.size1();j++)
					{
						boostgenotype_column(j,i) = _genotypes(j,i);
					}
				}
			}*/
			//We get the snp column number and get datas from it
			//boostgenotype_column = boostgenotype_column(current_SNP);
			G2_conditional_test_indep cond_g2(boostgenotype_column, _phenotypes, mb_temp);
			number_of_indep_test ++;

			// A G2 independency test is reliable when there are enough observations in each cell of the contingency table.
			// See the method Contingency::is_reliable() for details.
			if(cond_g2.is_reliable())
			{
				if(cond_g2.pval() < best_subset_pvalue)
				{
					best_subset_pvalue = cond_g2.pval();
					best_snp_index = i;
				}

				scores[i].push_back(cond_g2.g2());
			}
		}
	}
	// Append the best subset if alpha < threshold
	if(best_subset_pvalue < _params.aco_alpha)
	{
		vector<unsigned> best_subset = all_snps_combinations[best_snp_index];  //copy
		for(unsigned int i; i<best_subset.size(); i++)
		{
			mb.push_back(best_subset[i]);
		}

		for(vector<unsigned>::iterator it=best_subset.begin(); it != best_subset.end(); ++it)
		{
			snp_table.erase(it);
		}
	}
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
void Smmb_ACO::backward_phase(vector<unsigned int> &mb, vector<unsigned int> &snp_table)
{
	for (auto it =mb.begin(); it != mb.end(); it++)
	{
		unsigned int current_mb_elem = *it;
		vector<unsigned int> current_comb;
		unsigned int size = 3; //TODO maximum size of a combination, ask if we need to add it in parameters file
		vector<vector<unsigned int>> all_snps_combinations;
		vector<unsigned int> mb_minus_elem = mb;
		mb_minus_elem.erase(it);
		//We sort for the combination, not sure if it is usefull ...
		sort(mb_minus_elem.begin(),mb_minus_elem.end());
		//We run all combinations on the mb_minus_elem index
		Miscellaneous::combinator(mb_minus_elem,all_snps_combinations,size);
		Miscellaneous::link_comb_to_snp(mb_minus_elem,all_snps_combinations);

		for (unsigned int i=0; i<all_snps_combinations.size(); i++)
		{
			vector<unsigned int> current_combination = all_snps_combinations[i];
			blas_column boostgenotype_column(blas_column(_genotypes,current_mb_elem));
			G2_conditional_test_indep cond_g2(boostgenotype_column, _phenotypes, current_combination);
			number_of_indep_test ++;
			//If the pvalue of the test is above accepted alpha risk, discard this snp from MB
			if(cond_g2.pval() > _params.aco_alpha)
			{
				//If the SNP is removed from the MB, we won't test it over different combinations,
				//so we stop this loop instance and start the next one
				mb.erase(it);
				break;
			}

		}

		//New sampling phase but in the Markov blanket minus the current element


	}
}


void Smmb_ACO::best_mbs(vector<vector<unsigned int>> &mbs)
{
	  vector<vector<unsigned int>>::iterator uniq_it;
	  vector<vector<unsigned int>> uniq_mbs;
	  cout << "unique markov blankets:";

	  for (uniq_it=mbs.begin(); uniq_it!=mbs.end(); uniq_it++)
	  {
		  pair<map<vector<unsigned int>,unsigned int>::iterator,bool> current_mb;
		  current_mb = mbs_count.insert(pair<vector<unsigned int> ,unsigned int>(*uniq_it,0));
		  if (current_mb.second==false)
		  {
		    cout << "element "<< current_mb.first->second << " already existed";
		  }
	  }
	  for (auto it = mbs_count.begin(); it != mbs_count.end(); ++it )
	      {
		  for(unsigned int it_next = 0; it_next<mbs.size(); it_next++)
		 		  {
		 			  vector<unsigned int> mbs_temp = mbs[it_next];
		 			 if (it->first == mbs_temp)
		 			 {
		 				 it->second += 1;
		 			 }
		 		  }


	      }
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
