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
#include <omp.h>
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
	//Initialization of tau, eta and pdf values according to the parameters file.
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

/*
 *The run_ACO method call every function or void required to run the Smmb_ACO algorithm.
 */

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
		compute_cumulative_distrib_proba();

		//Parrallel computation of each ant. Have a possible memory leak that disappear when parrallel computing is removed.
//#pragma omp parallel for
		//Iterate over the number of ants in the parameters file
		for(unsigned int ant = 0; ant<_params.number_ants; ant++)
		{
			// Define a SNP set which is going to be selected by the sampling function
			vector<unsigned int> snp_table;
			//#pragma omp critical
			//Random picking of snps according to the ACO snp set size given in parameters file
			snp_sampling(snp_table);

			// learn MB from these SNPS
			list<unsigned int> mb;
			//learn mb is going to call forward and backward phase in order to learn a markov blanket.
			learn_mb(mb,snp_table);
			//Uncomment this step if you want the row list of learned markov blanket and call the best_mbs void.
			// If the markov blanket is not empty, add the learned markov blanket to the mbs vector.
			if(!mb.empty())
			{
				vector<unsigned int> mb_temp;
				//For every element in the mb list. We create a temporary
				for(auto i = mb.begin();i!= mb.end();i++)
				{
					mb_temp.push_back(*i);
				}
				mbs.push_back(mb_temp);
			}

		}
		//Update the pheromone levels.
		update_tau();

	}
}

/*
 *The sum_tau method compute the sum of tau with each value of the tau vector.
 */

void Smmb_ACO::sum_tau()
{
	// sum_of_tau initialisation at 0 (this void is called at each ACO iteration)
	sum_of_tau= 0.0;
	//Tried to use list type as requested for tau and eta variable but list don't use direct access, vector does.
	//We use each element of the tau vector to compute the sum of tau.
	for(unsigned int i=0; i<tau.size(); i++)
	{
		sum_of_tau += pow(tau.at(i), _params.aco_alpha) * pow(eta.at(i), _params.aco_beta);
	}

}

/*
 *The evaporation_rate_update method will update the tau value for a given snp according to its g2_score.
 */

void Smmb_ACO::evaporation_rate_update(unsigned int snp_index, float g2_score)
{
	tau[snp_index] = (1-_params.aco_rho) * tau[snp_index] + g2_score*_params.aco_lambda;
}

/*
 *The update_tau method is called at the end of an ACO iteration and will call the evaporation_rate_update method on each snp according to the
 *g2 score of this snp stored in the scores map.
 */


void Smmb_ACO::update_tau()
{
	//for each key in the scores map, we get the g2 score and use it in the evaporation_rate_update method
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

/*
 *The pheromone_for_snp method compute the effective probability for a given SNP and will be used to compute the distribution of probability for
 *an ACO iteration.
 */

float Smmb_ACO::pheromone_for_snp(float tau_for_snp, float eta_for_snp)
{
	float effective_proba = 0.0;
	tau_for_snp = pow(tau_for_snp, _params.aco_alpha)*pow(eta_for_snp, _params.aco_beta);

	effective_proba = tau_for_snp/sum_of_tau;

	return effective_proba;
}

/*
 *The compute_distrib_prob method call the pheromone_for_snp method for each tau score (linked to a SNP) in the tau vector.
 */

void Smmb_ACO::compute_distrib_prob()
{
	for(unsigned int i=0; i<tau.size(); i++)
	{
		pdf.at(i) = pheromone_for_snp(tau.at(i),eta.at(i));
	}
}

/*
 *The compute_cumulative_distrib_proba add the pdf value to the cumulated distribution of probabilities value if the pdf is greater than 0.
 */


void Smmb_ACO::compute_cumulative_distrib_proba()
{
	float memory = 0.0;
	cumulated_distrib_prob.clear();
	for(unsigned int i=0; i<pdf.size(); i++)
	{
		if (pdf.at(i) >0)
		{
			memory += pdf.at(i);
			//#pragma omp critical
			cumulated_distrib_prob[memory].push_back(i);
		}
	}
}

/*
 *The select_snp_in_distrib_prob method will select a SNP according to a threshold gave by a random proba.
 */

unsigned int Smmb_ACO::select_snp_in_distrib_prob(float prob)
{

	//it_after_prob tell the iterator to start after a key equal or greater than prob
	//it's like for (i=2; i<vector.size(); i++) we don't start at 0 but at lower_bound set at 2
	//map<float, vector<unsigned int>>::iterator it_after_prob;
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

/*
 *The snp sampling method is going to draw snps in the distribution of probabilites. It call the select_snp_on_distrib_prob.
 *If the snp is not already in the snp_table we add it to the selection
 */

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
			 //TODO : Memory leak seems to come from select_snp_in_distrib_prob. Struggling to debug this
			unsigned int drawn_snp = select_snp_in_distrib_prob(proba);
#pragma omp critical
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

/*
 *The learn_mb method is going to call the forward and backward in order to learn a markov blanket.
 */

void Smmb_ACO::learn_mb(list<unsigned int> &mb, vector<unsigned int> &snp_table)
{
	unsigned int counter = 0;
	list<unsigned int> memory_mb; //This a the vector of all the trials to learn a mb during this void
	//While the number of maximum trials to learn a mb is not reached or the mb is empty or not stabilized
	while((counter < _params.max_trials_learn_mb) && (mb.empty() && (memory_mb != mb)));
	{
		memory_mb = mb;
		forward_phase(mb,snp_table);
		backward_phase(mb,snp_table);
		counter++;
	}
	backward_phase(mb,snp_table);
}

/*
 *The forward_phase is going to statistically test the viability of the integration of a SNP in the markov blanket.
 */

void Smmb_ACO::forward_phase(list<unsigned int> &mb, vector<unsigned int> &snp_table)
{
	vector<unsigned int> random_snps;
	float best_subset_pvalue = 1.1;
	unsigned int best_snp_index = 100;
	unsigned int current_SNP = 0;
	//Sampling phase. SNPs are randomly drawn in the snp_table
	Miscellaneous::random_subset(snp_table,random_snps,_params.smallest_subset_size ,rand_seed);
	//Sorting is important for the combination method.
	sort(random_snps.begin(),random_snps.end());
	//Generating all combinations for the snp list
	vector<vector<unsigned int>> all_snps_combinations;
	//The combinator method search every combination in the drawn snps.
	Miscellaneous::combinator(random_snps,all_snps_combinations,_params.size);
	//This method link the previously made combination to the snp index.
	Miscellaneous::link_comb_to_snp(random_snps,all_snps_combinations);
	//for all combinations according to the size of the pattern, we pick a combination and
	//put it in the current combination vector.
	for (unsigned int x=0; x<all_snps_combinations.size(); x++)
	{
		vector<unsigned int> current_combination = all_snps_combinations.at(x);
		//Copy of the current combination vector into a temporary one because
		//we don't want to modify the vector being read in the following loop
		vector<unsigned int> current_combination_temp = current_combination;
		for(auto it =current_combination.begin(); it != current_combination.end(); it++)
		{
			current_SNP = *it;
			//We convert the markov blanket in a temporary vector type markov blanket.
			//The idea here is to use the direct acces properties of vectors in the next steps
			//This step of type conversion is very fast because markov blankets are small lists.
			vector<unsigned int> mb_temp;
			for(auto i = mb.begin();i!= mb.end();i++)
			{
				mb_temp.push_back(*i);
			}
			//erase current_SNP from temporary combination vector
			vector<unsigned int>::iterator INT_oui;
			unsigned int current_SNP_index;
			//The find method is going to seek the current SNP in the current combination_temp vector.
			INT_oui = find (current_combination_temp.begin(), current_combination_temp.end(), current_SNP);
			//If the current SNP is found, we erase it from the current_combination_temp vector
			if (INT_oui != current_combination_temp.end())
			{
				//Get the index of the SNP to erase
				current_SNP_index = distance(current_combination_temp.begin(), INT_oui);
				//erase goes from the beginning of the vector to the required indice.
				current_combination_temp.erase(current_combination_temp.begin()+current_SNP_index);
			}else
			{
				cout << "Element not found in current_combination_temp\n";
			}

			//Adding the current_combination without the current SNP to mb_temp
			for(unsigned int i=0; i<current_combination_temp.size(); i++)
			{
				mb_temp.push_back(current_combination_temp.at(i));
			}
			//Iterate over the current_combination_temp vector's size in order to test each
			//SNP of the combination against the MB currently being learned
			for(unsigned int i=0; i<current_combination_temp.size(); i++)
			{
				//Create a matrix of size number of patient and one SNP
				blas::matrix<int> boostgenotype_column(_genotypes.size1(),1);
				//The following loop iterate over the number of samples(patient) and get the current snp data in a blas matrix of one column
				for (unsigned int j = 0; j < _phenotypes.size1(); j++)
				{

					boostgenotype_column(j,0) =  _genotypes(j,current_SNP);
				}
				//If the mb is empty we go to the next iteration. We don't want to test a SNP conditionally to nothing.
				if(mb_temp.size()==0)
					break;
				//Run a G2 conditional independancy test with the currently evaluated SNP's data, phenotypes' data, the temporary mb
				//and the genotypes' data. Select true if you want to print contingency tables on stdout.
				G2_conditional_test_indep cond_g2(boostgenotype_column, _phenotypes, mb_temp,_genotypes, false);
				number_of_indep_test ++;

				// A G2 independency test is reliable when there are enough observations in each cell of the contingency table(more than 5).
				// See the method Contingency::is_reliable() for details.
				//If the pvalue of the above test is greater than the current best p_value for this subset, we had the current SNP and its score
				//to the scores map.
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
				//Iterate over the size of the vector containing the best_snp
				for(unsigned int j= 0; j<all_snps_combinations[best_snp_index].size();j++)
				{
					//add snps from the previously selected subset to the best_subset vector
					best_subset.push_back(all_snps_combinations.at(best_snp_index).at(j));
				}
				//If the best subset is not empty, we iterate over its size and add snps constituting the best subset to the mb
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

/*
 *The backward phase method is going to use a similar method to the forward phase but with the SNPs previously selected in forward phase
 *in order to check if snps in the current mb are still supposed to be a part of the markov blanket.
 */

void Smmb_ACO::backward_phase(list<unsigned int> &mb, vector<unsigned int> &snp_table)
{
	for(unsigned int i = 0; i<mb.size();i++)
	{
		//Iterator on the markov blanket
		std::list<unsigned int>::iterator it = mb.begin();
		//move forward the iterator with the i value.
		std::advance(it, i);
		//initialize the current_mb_elem. Take value of the iterator over the mb.
		unsigned int current_mb_elem = *it;
		vector<unsigned int> current_comb;
		vector<vector<unsigned int>> all_snps_combinations;
		vector<unsigned int> mb_minus_elem;
		//TODO: The list container is not the more optimized choice but is the working one.
		vector<unsigned int> mb_temp;
		//Copy markov blanket in a vector in order to use direct access
		for(auto i = mb.begin();i!= mb.end();i++)
		{
			mb_temp.push_back(*i);
		}
		//sort and erase duplicate values in the markov blanket.
		sort(mb_temp.begin(), mb_temp.end());
		mb_temp.erase(unique(mb_temp.begin(), mb_temp.end()), mb_temp.end());

		//add element to a temporary vector made of all the mb_temp except the element that we are going to evaluate with the G2_cond_test
		for(auto i = mb_temp.begin();i!= mb_temp.end();i++)
		{
			//If there is only one element in the markov blanket, we always add this element to the markov blanket and break the loop
			if(mb_temp.size()==1)
			{
				mb_minus_elem.push_back(*i);
				break;
				//if iterator value is not equal to the current_mb_elem that is going to be evaluated, add *i to mb_minus_elem
			}else if(*i != current_mb_elem)
			{
				mb_minus_elem.push_back(*i);
			}

		}
		//If the mb is empty, stop the backward phase.
		if(mb.size()==0)
		{
			break;
		}
		//We sort the mb_minus_elem before the combination phase.
		sort(mb_minus_elem.begin(),mb_minus_elem.end());
		//We run all combinations on the mb_minus_elem index
		Miscellaneous::combinator(mb_minus_elem,all_snps_combinations,_params.size);
		Miscellaneous::link_comb_to_snp(mb_minus_elem,all_snps_combinations);

		//Create a matrix of size number of patient and one SNP
		blas::matrix<int> boostgenotype_column(_genotypes.size1(),1);
		//Initialize a vector to store g2_results that will be add later to the results map
		vector<double> g2_results_temp;
		//Filling the matrix of genotypes datas for one SNP
		for (unsigned int j = 0; j < _phenotypes.size1(); j++)
		{

			boostgenotype_column(j,0) =  _genotypes(j,current_mb_elem);
		}
		// iterate over all the combinations made according to the size of epistasis and the number of snps.
		for (unsigned int i=0; i<all_snps_combinations.size(); i++)
		{
			vector<unsigned int> current_combination;
			//For element in a snp combination, we store them in a vector anmed current_combination.
			for(unsigned const& elem: all_snps_combinations[i])
			{
				current_combination.push_back(elem);
			}
			if(mb_temp.size()==0)
				break;
			//The following g2 conditional test of independancy will allow us to test if the SNP is relevant.
			G2_conditional_test_indep cond_g2(boostgenotype_column, _phenotypes, current_combination,_genotypes, false);
			number_of_indep_test ++;
			//We push pval and g2 score in the results vector
			g2_results_temp.push_back(cond_g2.pval());
			g2_results_temp.push_back(cond_g2.g2());
			//if the g2 is reliable, we add a one to the vector and if it's not we add a 0
			if(cond_g2.is_reliable())
			{
				g2_results_temp.push_back(1); //if reliable, add 1
			}else
			{
				g2_results_temp.push_back(0); // if no, add 0
			}
			//initialize the number of occurence
			g2_results_temp.push_back(1);

			map<vector<unsigned int>,vector<double>>::iterator it;
			//search for the current combination
			it = results.find(mb_temp);
			//if this combination exist in results, we add an occurence to the result.
			if (it != results.end())
			{
				#pragma omp critical //TODO to remove
				results[mb_temp][3] = results[mb_temp][3]+1;
			}else
			{
				//if this is the first occurency, we add the vector to the result with the proper key (mb_temp)
#pragma omp critical
				results[mb_temp] = g2_results_temp;
			}
			//If the g2 pvalue is above the threshold, we remove the SNP from the markov blanket and from the mb_temp.
			if(cond_g2.pval() > _params.alpha)
			{
				//If the SNP is removed from the MB, we won't test it over different combinations,
				//so we stop this loop instance and start the next one
				mb.remove(current_mb_elem);
				mb_temp.erase(remove(mb_temp.begin(), mb_temp.end(), current_mb_elem), mb_temp.end());
				//break; //TODO to remove
			}
		}
	}
}

/*
 *The best_mbs method is not usefull in the program but can be called in the main.cpp if you want to access to the list of markov blanket
 *that have been learnt during the Smmb_ACO iteration.
 */
void Smmb_ACO::best_mbs(vector<vector<unsigned int>> &mbs)
{
	//Iterate over the vector of mbs and then over the cotnent of a mb in order to print it in the stdout.
	for (unsigned int i = 0; i<mbs.size(); i++)
	{
		for(auto it = mbs_count[i].begin(); it != mbs_count[i].end(); it++)
		{
			cout << *it << " ";
		}
		cout << "##" << endl;
	}
	cout <<"#########################################################################"<<endl;
	cout << "fin de best mbs"<<endl;
}

/*
 *The write_results method write the results of Smmb_ACO in a file that is created where the output_path tells it to write.
 */

void Smmb_ACO::write_results(string output_path){
	//output file handling
	ofstream file;
	file.open(output_path);

	//sorting the results map according to the p_value. This step need a temporary vector called _optimum_set_vector.
	vector<pair<vector<unsigned>, vector<double>>> _optimum_set_vector;
	//iterate on the results map
	for (auto pattern : results)
	{
		//pair the key with its value and add it to the temporary vector called _optimum_set_vector
		pair<vector<unsigned>, vector<double>> pair_sort(pattern.first, pattern.second);
		_optimum_set_vector.push_back(pair_sort);
	}
	//sort the vector. Thanks to the pairing done in the previous loop, key and values are still linked. sorting according to p_values.
	//check the compareFunc for the details
	sort(_optimum_set_vector.begin(), _optimum_set_vector.end(), Miscellaneous::compareFunc);
	//if file is opened then write
	if(file)
	{
		//remove the following line if you want to remove the header.
		file << "Epistasis Pattern      p-value      score      reliable       occurences" <<endl;
		//iterate over the pairs in the optimum_set_vector.
		for(unsigned int i = 0 ; i<_optimum_set_vector.size() ; i++)
		{
			unsigned int j = 1;
			file << "{";
			//write the pattern in the file. we iterate on the value of the first element of the pair.
			for (auto key_it = _optimum_set_vector[i].first.cbegin(); key_it != _optimum_set_vector[i].first.cend(); key_it++)
			{
				// if it's not the last element of the pattern, add a comma after the SNP.
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
			//write the values linked to the first element of the pair of variables.
			for (auto val_it = _optimum_set_vector[i].second.cbegin(); val_it != _optimum_set_vector[i].second.cend(); val_it++)
			{
				//In order to improve the output, it could be useful to find a way to fix the columns in the print output.
				file << *val_it << "      ";
			}
			file << endl;
		}
		file.close();
	}else
		cerr << "Cannot open file." << endl;
}
