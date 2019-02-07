/*
 * SmmbACO.cpp
 *
 *  Created on: 9 nov. 2018
 *      Author: Louison Fresnais, François Courtin
 *      Project: SMMB-ACO and Genetic Algorithm for epistasis detection
 *      Under the supervision of Christine Sinoquet(Nantes University)
 *  Modified on: 05 fev 2018
 */

#include "SmmbACO.hpp"
using namespace std;

//Constructor
Smmb_ACO::Smmb_ACO(blas::matrix<int> & genos, blas::matrix<int> & phenos, Parameters_file_parsing params): _genotypes(genos), _phenotypes(phenos), _params(params)
{

	//Output file opening
    file_basename = basename((char*)_params.genos_file.c_str());

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

void Smmb_ACO::run_ACO()
{
	cout << "SMMB is currently running" << endl;
	for(unsigned number_aco_it=0; number_aco_it< _params.number_aco_iter; number_aco_it++)
	 {
		cout << "ACO iteration number: " << number_aco_it << endl;

		//Sum of pheromones
		sum_tau();
		//We clear the score map(key > value) for the new ACO iteration
		scores.clear();
		//Computation of the probability distribution
		compute_distrib_prob();
		compute_cumulative_dristrib_proba();

		//Parrallel computation of each ant

		#pragma omp parallel for
        for(unsigned int ant=0; ant<_params.number_ants; ant++)
        {
            // Define a SNP set which is going to be selected by the sampling function
            vector<unsigned int> snp_table;
            snp_sampling(snp_table);

            // learn MB from these SNPS
            list<unsigned int> mb;
            learn_mb(mb,snp_table);
            // add candidate mb to _mbs
            if(!mb.empty())
            {
    			vector<unsigned int> mb_temp;
    			for(auto i = mb.begin();i!= mb.end();i++)
    			{
    				mb_temp.push_back(*i);
    			}
            	mbs.push_back(mb_temp);
            }

        }
        update_tau(); //TODO

	 }
}


void Smmb_ACO::sum_tau()
{
    // sum_of_tau initialisation at 0 (this void is called at each ACO iteration)
	sum_of_tau= 0.0;
	//Tried to use list type for tau and eta variable but list don't use direct access, vector does.
    for(unsigned int i=0; i<tau.size(); i++)
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
		vector<double> scores_vec = it->second;
		for(unsigned int i = 0; i < scores_vec.size(); i++)
		{
			evaporation_rate_update(snp_index,scores_vec[i]);
		}
	}
}

float Smmb_ACO::pheromone_for_snp(float tau_for_snp, float eta_for_snp)
{
	float effective_proba = 0.0;
	tau_for_snp = pow(tau_for_snp, _params.aco_alpha)*pow(eta_for_snp, _params.aco_beta);

	effective_proba = tau_for_snp/sum_of_tau;

	return effective_proba;
}

//This method compute the distribution of probability
void Smmb_ACO::compute_distrib_prob()
{
	for(unsigned int i=0; i<tau.size(); i++)
	{
		pdf.at(i) = pheromone_for_snp(tau.at(i),eta.at(i));
	}
}

//This method ad to the cumulated distribution of probabilities table, every proba that is greater than 0
void Smmb_ACO::compute_cumulative_dristrib_proba()
{
	float memory = 0.0;
	cumulated_distrib_prob.clear();
	for(unsigned int i=0; i<pdf.size(); i++)
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
			float proba = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
			unsigned int drawn_snp = select_snp_in_distrib_prob(proba);

			//if nothing found, the find method return the last position of the vector
			if(find(snp_table.begin(), snp_table.end(), drawn_snp) == snp_table.end())
			{
				//if drawn_snp is not in snp_table, we add it
				snp_table.push_back(drawn_snp);
				snp_in_sample = true;

			}else
			{
				snp_in_sample = false;
			}
		}
		snp_in_sample = false;

	}

}

void Smmb_ACO::learn_mb(list<unsigned int> &mb, vector<unsigned int> &snp_table)
{
	unsigned int counter = 0;
	list<unsigned int> memory_mb; //This a the vector of all the trials to learn a mb during this void
	while((counter < _params.max_trials_learn_mb) && (mb.empty() && (memory_mb != mb)));
	{
		memory_mb = mb;
		forward_phase(mb,snp_table);
		backward_phase(mb,snp_table);
		counter++;
	}
	backward_phase(mb,snp_table);
}
void Smmb_ACO::forward_phase(list<unsigned int> &mb, vector<unsigned int> &snp_table)
{
	vector<unsigned int> random_snps;
	float best_subset_pvalue = 1.1;
	unsigned int best_snp_index = 100;
	unsigned int current_SNP = 0;
	//New sampling phase
	Miscellaneous::random_subset(snp_table,random_snps,_params.smallest_subset_size ,rand_seed);
	//Sorting is important for the combination method.
	sort(random_snps.begin(),random_snps.end());
	//Generating all combinations for the snp list
	vector<vector<unsigned int>> all_snps_combinations;
	Miscellaneous::combinator(random_snps,all_snps_combinations,_params.size);
	Miscellaneous::link_comb_to_snp(random_snps,all_snps_combinations);

	for (unsigned int x=0; x<all_snps_combinations.size(); x++)
	{
		vector<unsigned int> current_combination = all_snps_combinations.at(x);
	}
	for (unsigned int x=0; x<all_snps_combinations.size(); x++)
	{
		vector<unsigned int> current_combination = all_snps_combinations.at(x);
		vector<unsigned int> current_combination_temp = current_combination;
		for(auto it =current_combination.begin(); it != current_combination.end(); it++)
		{
			current_SNP = *it;
			//We merge the temporary markov blanket
			vector<unsigned int> mb_temp;
			for(auto i = mb.begin();i!= mb.end();i++)
			{
				mb_temp.push_back(*i);
			}
			//erase current_SNP from temporary combination vector
			vector<unsigned int>::iterator INT_oui;
			unsigned int current_SNP_index;
			INT_oui = find (current_combination_temp.begin(), current_combination_temp.end(), current_SNP);
			if (INT_oui != current_combination_temp.end())
			{
				current_SNP_index = distance(current_combination_temp.begin(), INT_oui);
				current_combination_temp.erase(current_combination_temp.begin()+current_SNP_index);
			}else
			{
				cout << "Element not found in current_combination_temp\n";
			}

			//Adding the current_combination whitout the current SNP to mb_temp
			for(unsigned int i=0; i<current_combination_temp.size(); i++)
			{
				mb_temp.push_back(current_combination_temp.at(i));
			}
			for(unsigned int i=0; i<current_combination_temp.size(); i++)
			{
				//Create a matrix of size number of patient and one SNP
				blas::matrix<int> boostgenotype_column(_genotypes.size1(),1);

				for (unsigned int j = 0; j < _phenotypes.size1(); j++)
				{

					boostgenotype_column(j,0) =  _genotypes(j,current_SNP);
				}
				if(mb_temp.size()==0)
					break;
				G2_conditional_test_indep cond_g2(boostgenotype_column, _phenotypes, mb_temp,_genotypes, false);
				number_of_indep_test ++;

				//We get the snp column number and get datas from it

				// A G2 independency test is reliable when there are enough observations in each cell of the contingency table.
				// See the method Contingency::is_reliable() for details.
				if(cond_g2.pval() < best_subset_pvalue)
				{
					vector<double> g2_results_temp;
					best_subset_pvalue = cond_g2.pval();
					best_snp_index = x;
					#pragma omp critical
					scores[current_SNP].push_back(cond_g2.g2());

				}
			}
		// Append the best subset if alpha < threshold
		if(best_subset_pvalue < _params.alpha)
		{
			vector<unsigned> best_subset;
			for(unsigned int j= 0; j<all_snps_combinations[best_snp_index].size();j++)
			{
				best_subset.push_back(all_snps_combinations.at(best_snp_index).at(j));
			}
			if(!best_subset.empty())
			{
				for(unsigned int i=0; i<best_subset.size(); i++)
				{
					list<unsigned int>::iterator findIter = find(mb.begin(), mb.end(), best_subset.at(i));
					if(findIter == mb.end())
					{
						mb.push_back(best_subset.at(i)); // Add snps from best_subset to the MB is not already find
					}
					//Erase best_subset values from the random_snps vector
					random_snps.erase(remove(random_snps.begin(), random_snps.end(), best_subset.at(i)), random_snps.end());
				}

			}

		}

	}
}
}
void Smmb_ACO::backward_phase(list<unsigned int> &mb, vector<unsigned int> &snp_table)
{

	for(unsigned int i = 0; i<mb.size();i++)
	{

		//TESTING ALTERNATIVE METHOD
		std::list<unsigned int>::iterator it = mb.begin();
		std::advance(it, i);
		unsigned int current_mb_elem = *it;
		//END OF ALTERNATIVE METHOD
		vector<unsigned int> current_comb;
		vector<vector<unsigned int>> all_snps_combinations;
        vector<unsigned int> mb_minus_elem;
    	//TODO: The list container is not the more optimized choice but is the working one atm
		vector<unsigned int> mb_temp;
		for(auto i = mb.begin();i!= mb.end();i++)
		{
			mb_temp.push_back(*i);
		}
		//Erase duplicate values in the markov_blanket
        sort(mb_temp.begin(), mb_temp.end());
        mb_temp.erase(unique(mb_temp.begin(), mb_temp.end()), mb_temp.end());


        for(auto i = mb_temp.begin();i!= mb_temp.end();i++)
        {
        	if(mb_temp.size()==1)
        	{
        		mb_minus_elem.push_back(*i);
        		break;
        	}else if(*i != current_mb_elem)
        	{
        		mb_minus_elem.push_back(*i);
        	}

        }
        if(mb.size()==0)
        {
        	break;
        }
		sort(mb_minus_elem.begin(),mb_minus_elem.end());
		//We run all combinations on the mb_minus_elem index
		Miscellaneous::combinator(mb_minus_elem,all_snps_combinations,_params.size);
		Miscellaneous::link_comb_to_snp(mb_minus_elem,all_snps_combinations);

        //Create a matrix of size number of patient and one SNP
        blas::matrix<int> boostgenotype_column(_genotypes.size1(),1);
        vector<double> g2_results_temp;
        for (unsigned int j = 0; j < _phenotypes.size1(); j++)
        {

        	boostgenotype_column(j,0) =  _genotypes(j,current_mb_elem);
        }
        for (unsigned int i=0; i<all_snps_combinations.size(); i++)
        {
        	vector<unsigned int> current_combination;
        	for(unsigned const& elem: all_snps_combinations[i])
        	{
        		current_combination.push_back(elem);
        	}
        	if(mb_temp.size()==0)
        		break;
        	G2_conditional_test_indep cond_g2(boostgenotype_column, _phenotypes, current_combination,_genotypes, false);
        	number_of_indep_test ++;
        	//If the pvalue of the test is above accepted alpha risk, discard this snp from MB
        	//TODO TRYING TO ADD THE RESULTS IN THE MARKOV BLANKET
        	g2_results_temp.push_back(cond_g2.pval());
        	g2_results_temp.push_back(cond_g2.g2());
        	if(cond_g2.is_reliable())
        	{
        		g2_results_temp.push_back(1); //if reliable, add 1
        	}else
        	{
        		g2_results_temp.push_back(0); // if no, add 0
        	}
        	g2_results_temp.push_back(0); // init number of occurences at 0 (not sure if it's useful
        	map<vector<unsigned int>,vector<double>>::iterator it;
        	it = results.find(mb_temp); //search for the current combination
        	if (it != results.end()) //if this combination exist in results
        	{
        		results[mb_temp][3] = results[mb_temp][3]+1;
        	}else
        	{
				#pragma omp critical
        		results[mb_temp] = g2_results_temp;
        	}
        	if(cond_g2.pval() > _params.alpha)
        	{
        		//If the SNP is removed from the MB, we won't test it over different combinations,
        		//so we stop this loop instance and start the next one
        		mb.remove(current_mb_elem);
        		mb_temp.erase(remove(mb_temp.begin(), mb_temp.end(), current_mb_elem), mb_temp.end());
        		break;
        	}

        }

        //New sampling phase but in the Markov blanket minus the current element


	}
}


void Smmb_ACO::best_mbs(vector<vector<unsigned int>> &mbs)
{
		vector<vector<unsigned int>>::iterator uniq_it;
		for (uniq_it=mbs.begin(); uniq_it!=mbs.end(); uniq_it++)
		{
			vector<unsigned int>temp_mb = *uniq_it;
			mbs_count.push_back(temp_mb);
		}
		for (unsigned int i = 0; i<mbs_count.size(); i++)
		{
			for(auto it = mbs_count[i].begin(); it != mbs_count[i].end(); it++)
			{
				cout << *it << " ";
			}
			cout << "##" << endl;
		}
		cout <<"#########################################################################"<<endl;
		cout <<"#########################################################################"<<endl;
		cout << "fin de best mbs"<<endl;
}
void Smmb_ACO::write_results(string output_path){
	ofstream file;
	file.open(output_path);

	//file << "number of markov blankets learnt: " << mbs.size()<<endl;
	//file << "number of different markov blankets in the results map" << results.size()<<endl;
	/* Map order*/
	vector<pair<vector<unsigned>, vector<double>>> _optimum_set_vector;
	for (auto pattern : results)
	{
	    pair<vector<unsigned>, vector<double>> pair_sort(pattern.first, pattern.second);
	    _optimum_set_vector.push_back(pair_sort);
	}
	sort(_optimum_set_vector.begin(), _optimum_set_vector.end(), Miscellaneous::compareFunc);

	if(file)
	{
		file << "Epistasis Pattern      p-value      score      reliable       occurences" <<endl;
		//TODO : Link the SNP index with the Phenotype header !
		for(unsigned int i = 0 ; i<_optimum_set_vector.size() ; i++)
		{
			unsigned int j = 1;
			file << "{";
			for (auto key_it = _optimum_set_vector[i].first.cbegin(); key_it != _optimum_set_vector[i].first.cend(); key_it++)
			{
				if(j < _params.size)
				{
					file << *key_it << ",";
				}else
				{
					file << *key_it;
				}

				j++;
			}
			file << "}        ";
			for (auto val_it = _optimum_set_vector[i].second.cbegin(); val_it != _optimum_set_vector[i].second.cend(); val_it++)
			{
				file << *val_it << "      ";
			}
			file << endl;
		}

		file.close();
	}else
		cerr << "Cannot open file." << endl;
}
