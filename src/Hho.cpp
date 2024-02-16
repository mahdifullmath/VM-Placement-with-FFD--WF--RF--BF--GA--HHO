/*
 * Hho.cpp
 *
 *  Created on: Jan 15, 2024
 *      Author: mahdi
 */
#include <cmath>
#include <time.h>
#include <algorithm>
#include <random>
#include <iostream>
#include <vector>
#include "Hho.h"

Hho::Hho()
{
}
Hho::~Hho()
{
}

Hho::Hho(int epoch, int pop_size, int num_dims)
{
	this->epoch = epoch;
	this->pop_size = pop_size;
	this->sort_flag = true;
	this->dimPlacement = num_dims;
	this->energy = max_energy_Hho;
	//this->Hhoin(pop_size);
}
/*void Hho::Hhoin( int pop_size) {

        pop.resize(pop_size);
        for(auto& solution : pop) {
            solution = Placement();
            random_shuffle(solution.structure.begin(), solution.structure.end());
        }
//        bestSolution = pop[0];
    }*/

void Hho:: run( vector<vm> vms,  vector<pm> servers) {

            //scoutPhase(vms, servers);
            if (energy > 1)			//exploration
            	explorePhase(vms, servers);
            else if (energy > 0.5)	//exploitation Soft besiege with progressive rapid dives
            	exploitPhase(vms, servers);
            else 					//exploitation Hard besiege with progressive rapid dives
            	attackPhase(vms, servers);

    }
/*
void Hho::scoutPhase(vector<vm> vms, vector<pm> servers) {
    for (int i = 0; i < pop_size; i++) {
        pop[i].calculateFitPlacement(servers, vms);
        if (pop[i].fit < bestFit) {
            bestFit = pop[i].fit;
            bestPlacement = pop[i];
        }
    }
}
 */
void Hho:: scoutPhase( vector<vm> vms,  vector<pm> servers) {
    for (int i = 0; i < pop_size; i++) {
        pop[i].insertVMListRandomHho(pop[i].findMissingVMIndicesHho(pop[i]));
        pop[i].checkValidityPlacement();
    }
 }
/*
void Hho::explorePhase(vector<vm> vms, vector<pm> servers) {
    for (int i = 0; i < pop_size; i++) {
    	//pop[i].mutate(pop[i]); // Assuming IndicesWithDifferentCPUHho is a member function of the Hho class
    	pop[i].explore();
        pop[i].updateFitness(); // Assuming there is a function to update the fitness of the placement
    }
    // Sort the population based on fitness
    this->sortPopFitnessDecreasing();

    this->setBestAndWorstPlacements(); // Assuming there is a function to update the best placement
}*/
/*
void Placement::explore() {
    // Perform the explore phase logic here
    // For example, you can implement a simple exploration strategy
    // such as randomly changing the placement of a VM within the structure

    // Select a random index within the structure
    int randIndex = rand() % structure.size();

    // Randomly select a new position for the VM within the structure
    int randNewIndex = rand() % structure.size();

    // Move the selected VM to the new position within the structure
    vmPlacement temp = structure[randIndex];
    structure[randIndex] = structure[randNewIndex];
    structure[randNewIndex] = temp;

    // Recalculate the fitness of the placement after the VM move
    calculateFitPlacement(PMList, VMList);

    // Update the best and worst fitness if necessary
    // Assuming bestFit and bestPlacement are member variables of the Placement class
    if (fit < bestFit) {
        bestFit = fit;
        bestPlacement = *this;  // Assuming bestPlacement is of type Placement
    }
}*/

void Hho::explorePhase(vector<vm> vms, vector<pm> servers) {
    // Iterate through the population
    for (int i = 0; i < pop_size - 1 ; i+=2) {
        // Select a random index within the structure
    	Placement newPlacement;
    	newPlacement = explore(pop[i] , pop[i+1] , servers , vms);

		int rand1 = rand() % (num_vm_experiment);
		int rand2 = rand() % (num_vm_experiment);
		if (rand2 < rand1)
			swap(rand2 , rand1);
		for (int j = rand1 ; j <= rand2 / 2 ; j++){
			swap(pop[i+1].structure[j], pop[i+1].structure[rand2-j]);
		}
    	// Recalculate the fitness of the placement after the VM move
    	pop[i+1].calculateFitPlacement(servers, vms);
    	updateBestPlacement(pop[i+1]);
    	updateWorstPlacement(pop[i+1]);
    	pop[i] = newPlacement;
        // Recalculate the fitness of the placement after the VM move
        pop[i].calculateFitPlacement(servers, vms);
        updateBestPlacement(pop[i]);
        updateWorstPlacement(pop[i]);
        // Update the best and worst fitness if necessary


    }
    //setBestAndWorstPlacements();
}

/*
void Hho::explorePhase( vector<vm> vms,  vector<pm> servers) {
        for(auto& solution : pop) {
            double fit = solution.fit;

            for(int i = 0; i < solution.structure.size(); ++i) {
                vector<int> neighbors;

                for(int j = 0; j < servers.size(); ++j) {
                    if(j != solution.structure[i].id) {
                        neighbors.push_back(j);
                    }
                }

                if(!neighbors.empty()) {
                    int index = generator() % neighbors.size();
                    int neighbor = neighbors[index];
                    solution.structure[i].id = neighbor;
                }
            }
        }
    }
    */
/*
void Hho::exploitPhase(vector<vm> vms, vector<pm> servers) {
    for (int i = 0; i < pop_size; i++) {
        pop[i].exploit(vms, servers); // Assuming there is an exploit function in the Placement class
        pop[i].updateFitness(); // Assuming there is a function to update the fitness of the placement
    }
    // Sort the population based on fitness
    this->sortPopFitnessDecreasing();

    this->setBestAndWorstPlacements(); // Assuming there is a function to update the best placement
}
void Hho::exploitPhase(vector<vm> vms, vector<pm> servers) {
    // Select a random placement from the population
    int randIndex = rand() % pop_size;

    // Calculate the fitness of the selected placement
    pop[randIndex].calculateFitPlacement(servers, vms);

    // Update the best and worst fitness if necessary
    if (pop[randIndex].fit < bestFit) {
        bestFit = pop[randIndex].fit;
        bestPlacement = pop[randIndex];
    }

    // Perform exploit by moving VMs within the placement to improve fitness
    moveVmsToImproveFitness(pop[randIndex], servers, vms);
}
*/
/*
void Hho::moveVmsToImproveFitness(Placement &p, vector<pm> servers, vector<vm> vms) {
    // Calculate the current fitness of the placement
    double currentFitness = p.calculateFitPlacement(servers, vms);

    // Iterate through VMs in the placement and attempt to improve fitness
    for (int i = 0; i < p.num_vms; i++) {
        // Consider moving the current VM to another server within the placement
        for (int j = 0; j < p.num_pms; j++) {
            if (j != p.vms[i].pm_id) {
                // Temporarily move the VM to the new server
                int originalPmId = p.vms[i].pm_id;
                p.vms[i].pm_id = j;

                // Recalculate fitness after the VM move
                double newFitness = p.calculateFitPlacement(servers, vms);

                // If the fitness improves, keep the VM on the new server
                if (newFitness < currentFitness) {
                    currentFitness = newFitness;
                    break;  // Move to the next VM
                } else {
                    // Revert the VM back to its original server
                    p.vms[i].pm_id = originalPmId;
                }
            }
        }
    }
}

*/

void Hho::exploitPhase( vector<vm> vms,  vector<pm> servers) {
    // Iterate through the population
    for (int i = 0; i < pop_size; i++) {
        // Select a random index within the structure
    	srand(time(0));
    	double random_number = (double)rand() / RAND_MAX;
    	if (random_number < Hho_expolit_rate){
    		//Placement newPlacement;
    		//newPlacement = explore(pop[i] , bestPlacement , servers , vms);
    		Placement newPlacement2 = pop[i];
    		int rand1 = rand() % (num_vm_experiment);
    		int rand2 = rand() % (num_vm_experiment);
    		swap(newPlacement2.structure[rand1], newPlacement2.structure[rand2]);
    		//pop[i] = newPlacement;
    		pop[i] = newPlacement2;
    		// Recalculate the fitness of the placement after the VM move
    		pop[i].calculateFitPlacement(servers, vms);
    		updateBestPlacement(pop[i]);
    		updateWorstPlacement(pop[i]);

    	}
    	else{
    		//pop[i] = bestPlacement;
    		int rand1 = rand() % (num_vm_experiment);
    		int rand2 = rand() % (num_vm_experiment);
    		if (rand2 < rand1)
    			swap(rand2 , rand1);
    		for (int j = rand1 ; j <= rand2 / 2 ; j++){
    			swap(pop[i].structure[j], pop[i].structure[rand2-j]);
    		}
    	}
		// Recalculate the fitness of the placement after the VM move
		pop[i].calculateFitPlacement(servers, vms);
		updateBestPlacement(pop[i]);
		updateWorstPlacement(pop[i]);

    }
    //setBestAndWorstPlacements();
}

void Hho::attackPhase( vector<vm> vms,  vector<pm> servers) {
	for (int i = 0; i < pop_size; i++) {

		Placement randPlacement = pop[i];
		Placement randPlacement2 = pop[i];
		int rand1 = rand() % (num_vm_experiment);
		int rand2 = rand() % (num_vm_experiment);
		swap(randPlacement2.structure[rand1], randPlacement2.structure[rand2]);
		if (rand2 < rand1)
		    swap(rand2 , rand1);
		for (int j = rand1 ; j <= rand2 / 2 ; j++){
		    swap(randPlacement.structure[j], randPlacement.structure[rand2-j]);
		    		}
//		swap(randPlacement.structure[rand1], randPlacement.structure[rand2]);
		randPlacement2.calculateFitPlacement(servers, vms);
		randPlacement.calculateFitPlacement(servers, vms);
		if(randPlacement.fit > randPlacement2.fit)
			swap(randPlacement,randPlacement2);
		float random_number = (float)rand() / RAND_MAX;
		if (randPlacement.fit < bestFit || random_number <  Hho_attack_rate){
			pop[i] = randPlacement;

		}
		else{
			pop[i] = bestPlacement;
			int rand1 = rand() % (num_vm_experiment);
			int rand2 = rand() % (num_vm_experiment);
			swap(randPlacement.structure[rand1], randPlacement.structure[rand2]);
		}
		pop[i].calculateFitPlacement(servers, vms);
		updateBestPlacement(pop[i]);
		updateWorstPlacement(pop[i]);
	}
    //setBestAndWorstPlacements();
}
/*
void Hho::attackPhase(vector<vm> vms, vector<pm> servers) {
    // Select two random placements from the population
    int randIndex1 = rand() % pop_size;
    int randIndex2 = rand() % pop_size;

    // Check if the two selected placements are different
    while (randIndex2 == randIndex1) {
        randIndex2 = rand() % pop_size;
    }

    // Perform attack by putting bubble in indices with different CPU
    putBubbleInIndicesWithDifferentCPUHho(pop[randIndex1], pop[randIndex2]);

    // Calculate fitness of the modified placements
    pop[randIndex1].calculateFitPlacement(servers, vms);
    pop[randIndex2].calculateFitPlacement(servers, vms);

    // Update the best and worst fitness if necessary
    if (pop[randIndex1].fit < bestFit) {
        bestFit = pop[randIndex1].fit;
        bestPlacement = pop[randIndex1];
    }
    if (pop[randIndex2].fit < bestFit) {
        bestFit = pop[randIndex2].fit;
        bestPlacement = pop[randIndex2];
    }
}
*/
/*
void Hho::putBubbleInIndicesWithDifferentCPUHho(Placement &p1, Placement &p2) {
    // Calculate the average CPU usage of the two placements
    double avgCpuUsage = (p1.totalCpuUsage() + p2.totalCpuUsage()) / 2.0;

    // Find VMs in p1 with CPU usage greater than the average
    vector<int> highCpuVmsP1;
    for (int i = 0; i < p1.num_vms; i++) {
        if (p1.vms[i].cpu_usage > avgCpuUsage) {
            highCpuVmsP1.push_back(i);
        }
    }

    // Find VMs in p2 with CPU usage less than the average
    vector<int> lowCpuVmsP2;
    for (int i = 0; i < p2.num_vms; i++) {
        if (p2.vms[i].cpu_usage < avgCpuUsage) {
            lowCpuVmsP2.push_back(i);
        }
    }

    // Perform VM swaps between p1 and p2 to balance CPU usage
    for (int i = 0; i < min(highCpuVmsP1.size(), lowCpuVmsP2.size()); i++) {
        // Swap the VMs between p1 and p2
        VM temp = p1.vms[highCpuVmsP1[i]];
        p1.vms[highCpuVmsP1[i]] = p2.vms[lowCpuVmsP2[i]];
        p2.vms[lowCpuVmsP2[i]] = temp;
    }
}
*/

//
//int randomNM(int start, int end)
//{
//	srand(time_t);
//	int range = (end - start) + 1;
//	int random_int = start + (rand() % range);
//	return random_int;
//}
//
void Hho::initialize(vector<vm> VMList)
{
	// create Placements where each Placement has a structure to hold the VMP solution
	for (unsigned int i = 0; i < pop_size; i++)
	{
		Placement Placement1;
		unsigned int rem_dim_Placement = dimPlacement;
		std::vector<unsigned int> vm_list;
		for (unsigned int i = 0; i < dimPlacement; i++)
		{
			vm_list.push_back(i);
		}
		// std::cout<<"initializing the pop"<<std::endl;
		for (unsigned int j = 0; j < dimPlacement; j++)
		{
			unsigned int Irand = rand() % vm_list.size();
			vmPlacement myVM;
			myVM.id = vm_list[Irand];
			myVM.vmInfo = VMList[Irand];
			Placement1.structure.push_back(myVM);
			// pop_temp[i].push_back(vm_list[Irand]);
			vm_list.erase(vm_list.begin() + Irand);
			Placement1.PlacementID = i;
		}
		pop.push_back(Placement1);
	}
	bestFit = 1e+20;
	worstFit = 0;
	bestStabLevel = 1e+30;
	bestPow_dis = 1e+30;
	bestResourceWastage = 1e+30;
	currentIteration = 0;
	crossOverCount = 0, mutationCount = 0;
	numFuncEval_Hho = 0;
}

void Hho::printPop()
{
	// std:: vector< std::vector<unsigned int> >::iterator row;
	// std:: vector<unsigned int>::iterator col;
	int i = 0;
	for (std::vector<Placement>::iterator it1 = pop.begin(); it1 != pop.end(); ++it1)
	{
		std::cout << "Placement-" << it1->PlacementID << ":";
		for (std::vector<vmPlacement>::iterator it2 = it1->structure.begin(); it2 != it1->structure.end(); ++it2)
		{
			std::cout << " " << it2->id;
		}
		i++;
		std::cout << std::endl;
	}
}
// TODO: assign a proper and valid ID to each newly created Placement
void Hho::printPlacementVector(vector<Placement> vectPlacements)
{
	cout << endl;
	for (std::vector<Placement>::iterator it1 = vectPlacements.begin(); it1 != vectPlacements.end(); ++it1)
	{
		std::cout << "Placement-" << it1->PlacementID << ":";
		for (std::vector<vmPlacement>::iterator it2 = it1->structure.begin(); it2 != it1->structure.end(); ++it2)
		{
			std::cout << " " << it2->id;
		}
		std::cout << std::endl;
	}
}

void Hho::calculateFitnessPop(vector<pm> PMList, vector<vm> VMList)
{
	int i = 0;
	for (auto xx : pop)
	{
		xx.calculateFitPlacement(PMList, VMList);
		if (xx.fit < bestFit)
		{
			bestPlacement = xx;
			bestFit = xx.fit;
		}
		if (xx.fit > worstFit)
		{
			worstPlacement = xx;
			worstFit = xx.fit;
		}
		if (VERBOSE)
			cout << "Atm-" << xx.PlacementID << " is " << xx.fit << endl;
		this->pop[i].fit = xx.fit;
		i++;
	}
	if (VERBOSE)
		cout << endl;
}

void Hho::sortPopFitnessDecreasing()
{
	// NOTE: a higher fitness is worse since power is used as the fitness.
	int i = 0;
	for (unsigned int i = 0; i < pop.size(); i++)
	{
		for (unsigned int j = i + 1; j < pop.size(); j++)
		{
			unsigned int tmp1 = pop[i].fit;
			unsigned int tmp2 = pop[j].fit;

			if (tmp1 > tmp2)
			{
				swap(pop[i], pop[j]);
			}
		}
	}
}

void Hho::printFitPop()
{
	for (std::vector<Placement>::iterator it1 = pop.begin(); it1 != pop.end(); ++it1)
	{
		std::cout << "Fit Placement-" << it1->PlacementID << ": " << it1->fit << endl;
	}
}

void Hho::evolve(vector<pm> PMList, vector<vm> VMList, std::ofstream &fileName, unsigned int iterationForAvg)
{
	char name_of[100];
	ofstream convergence_file_Hho;
	string str_of1;
	string nvmstr = std::to_string(num_vm_experiment);
	string pselstr = std::to_string(EXPERIMENT_P_VALUE_SELECTOR);
	sprintf(name_of, "%s/CONV_%s_%sVMs_%sPvalue.txt", "results", "Hho", nvmstr.c_str(), pselstr.c_str());
	str_of1 = name_of;

	if (iterationForAvg == 0)
	{
		convergence_file_Hho.open(name_of, std::ofstream::out | std::ofstream::trunc);
		convergence_file_Hho.close();
	}
	convergence_file_Hho.open(str_of1, fstream::app);
	std::vector<std::vector<double>> pop_new;
	vector<Placement> newPop;
	unsigned int cnt_print = 0;
	currentIteration = 0;
	while (currentIteration < Hho_epoch and numFuncEval_Hho <= max_FE_Hho)
	{ // unsigned int idx = 0; idx < this->epoch; idx++) {
		// calculate fitness of all Placements
		// for(std::vector<Placement>::iterator it1 = pop.begin(); it1 != pop.end(); ++it1){
		newPop.erase(newPop.begin(), newPop.end()); // newPop is formed during each iteration
		// cout<<"best / worst fits"<<bestFit<<" / "<<worstFit<<endl;
		this->calculateFitnessPop(PMList, VMList);
		this->sortPopFitnessDecreasing();
		if (VERBOSE)
			this->printFitPop();
		this->modifyPlacementIndices();
		// cout<<"best / worst fits"<<bestFit<<" / "<<worstFit<<endl;
		// this->setBestAndWorstPlacements();
		// cout<<"best / worst fits"<<bestFit<<" / "<<worstFit<<endl;
		convergence_file_Hho << currentIteration << "  " << bestFit << endl;
		if (VERBOSE)
			printPop();
		// cout<<"Best fit in iteration-"<<idx<<": "<<bestFit<<endl;

//		int s = (10 * pop_size) / 100;
//		for (int i = 0; i < s; i++)
//			newPop.push_back(pop[i]);
//
//		// From 50% of fittest pop, Placements
		// will mate to produce offspring
//		s = pop_size - s; //(90*pop_size)/100;
//		for (int i = 0; i < s; i++)
//		{
//			//   int len = pop.size();
//			int r = randomNM(0, pop.size() - 1);
//			Placement parent1 = pop[r];
//			r = randomNM(0, pop.size() - 1);
//			Placement parent2 = pop[r];
//			//Placement offspring = CrossOver(parent1, parent2, PMList, VMList);
//			//newPop.push_back(offspring);
//		}
//		newPop = mutate(newPop);
		run(VMList, PMList);
		currentIteration++;
		energy = energy - max_energy_Hho / Hho_epoch;
		cnt_print++;
		// Hho_convergence << cnt<<","<<pop[0].fit<<"\n";
		cout <<"(HHO):" <<currentIteration << "...";
		if (cnt_print == 100)
		{
			cout << endl;
			cnt_print = 0;
		}

		if (VERBOSE)
			printPlacementVector(addedPlacementsIteration);
//		pop.clear();
		// Placements added/removed in alpha/beta/Hhomma reaction are now added/removed to pop
		// note: in Positron, the worst Placements are merged and added to pop while the worst Placements are removed
		if (VERBOSE)
			printPlacementVector(addedPlacementsIteration);
//		for (auto z : addedPlacementsIteration)
//			pop.push_back(z);
		addedPlacementsIteration.clear(); //(addedPlacementsIteration.begin(),addedPlacementsIteration.end());
		/* for(auto z:toBeRemovedPlacementsIteration){
			unsigned int idxPlacement=findPlacementIDxInPop(z);
			pop.erase(pop.begin()+idxPlacement);
		}*/
		this->modifyPlacementIndices();
		toBeRemovedPlacementsIteration.clear(); //(toBeRemovedPlacementsIteration.begin(),toBeRemovedPlacementsIteration.end());
		cout << "Pop size:" << pop.size() << endl;
	}
	float power = bestPlacement.calculatePowerPlacement(PMList, VMList);

	Solution = bestPlacement.fitPlacementSequenceUsingFF(PMList, VMList);
	convergence_file_Hho.close();
	// Solution.numActivePMs=num_act_pms;
	Solution.power = power;
	Solution.numVMs = VMList.size();
	// Solution.calculateMappingPower(PMList);
	Solution.calculateMappingWastageAndBalance(PMList);
	// cout<<"Out LB CPU: "<<Solution.loadBalanceMEM<<endl;
	//  cout<<"Out LB MEM: "<<Solution.loadBalanceCPU<<endl;
	//  Solution.calculateMappingLodaBalance(PMList);
	Solution.iterForAvg = iterationForAvg;
	Solution.printMappingStatsIntoFile(fileName, PMList);
}

vector<Placement> Hho::sortPopByFitnessDecreasingly()
{
 /*   vector<Placement> sortedPop = pop;
    sort(sortedPop.begin(), sortedPop.end(),
        [](const Placement& a, const Placement& b) {
            return a.fit < b.fit;
        });*/
	vector<Placement> sortedPop;
	sortedPop = pop;
	for (unsigned int i = 0; i < sortedPop.size(); i++)
	{
		for (unsigned int j = i + 1; j < sortedPop.size(); j++)
		{
			float temp1 = sortedPop[i].fit;
			float temp2 = sortedPop[j].fit;
			if (temp1 < temp2)
			{
				swap(sortedPop[i], sortedPop[j]);
			}
		}
	}
	return sortedPop;
}

unsigned int Hho::findPlacementIdInPop(Placement thePlacement)
{
	unsigned int idx = 0;
	for (auto x : pop)
	{
		if (x.PlacementID == thePlacement.PlacementID)
		{
			return idx;
		}
		idx++;
	}
	cout << endl
		 << "ATOM INDEX(" << thePlacement.PlacementID << ") NOT FOUND!!!, CHECK IT" << endl;
	return 0;
}

unsigned int Hho::findValidPlacementID()
{
	vector<unsigned int> allPlacementIndices(pop.size() + 1000);
	vector<unsigned int> missingPlacementID;
	iota(allPlacementIndices.begin(), allPlacementIndices.end(), 0);
	for (auto x : allPlacementIndices)
	{
		bool exists = false;
		for (auto y : pop)
		{
			if (y.PlacementID == x)
			{
				exists = true;
				break;
			}
		}
		for (auto y : addedPlacementsIteration)
		{
			if (y.PlacementID == x)
			{
				exists = true;
				break;
			}
		}
		if (!exists)
		{
			return x;
		}
	}
	cout << "ERRO: no proper Placement ID found!!";
	return 0;
}

void Hho::updateBestPlacement(Placement theNewPlacement)
{
	if (bestFit > theNewPlacement.fit)
	{
		bestPlacement = theNewPlacement;
		bestFit = theNewPlacement.fit;
	}
}

void Hho::updateWorstPlacement(Placement theNewPlacement)
{
	if (worstFit < theNewPlacement.fit)
	{
		worstPlacement = theNewPlacement;
		worstFit = theNewPlacement.fit;
	}
}

void Hho::modifyPlacementIndices()
{
	unsigned int PlacementID = 0;
	for (auto x : pop)
	{
		pop[PlacementID].PlacementID = PlacementID;
		PlacementID++;
	}
}

void Hho::setBestAndWorstPlacements()
{
	float fitArray[pop.size()];
	unsigned int fitIndexInPop[pop.size()];
	unsigned int i = 0;
	for (auto x : pop)
	{
		// updateBestPlacement(x);
		// updateWorstPlacement(x);
		fitArray[i] = x.fit;
		fitIndexInPop[i] = i; // hold i-th index in Pop instead of PlacementID
		i++;
	}

	for (unsigned int ii = 0; ii < pop.size(); ii++)
	{
		for (unsigned int jj = 0; jj < pop.size(); jj++)
		{
			if (fitArray[ii] < fitArray[jj])
			{
				float fitTemp = fitArray[ii];
				fitArray[ii] = fitArray[jj];
				fitArray[jj] = fitTemp;
				fitIndexInPop[ii] = fitIndexInPop[jj];
			}
		}
	}
	/* cout<<"Fit values:"<<endl;
	 for(unsigned int ii=0;ii<pop.size();ii++){
		 cout<<" "<<fitArray[ii];
	 }
	 cout<<endl;*/
	worstPlacement = pop.at(fitIndexInPop[pop.size() - 1]);
	worstFit = fitArray[pop.size() - 1];
	bestPlacement = pop.at(fitIndexInPop[0]);
	bestFit = fitArray[0];
}

vector<Placement> Hho::mutate(vector<Placement> popPlacements)
{
	mutationCount++;
	int mut_num = (int)MUTATION_RATE * dimPlacement * pop_size;
	for (int i = 0; i < mut_num; i++)
	{
		int rnd = rand() % popPlacements.size();
		Placement randPlacement = popPlacements[rnd];
		int rand1 = rand() % (num_vm_experiment);
		int rand2 = rand() % (num_vm_experiment);
		swap(popPlacements[rnd].structure[rand1], popPlacements[rnd].structure[rand2]); // child_Placement.push_back(mutated_Placements());
	}
	return popPlacements;
}
Placement Hho::explore(Placement par1, Placement par2, vector<pm> PMList, vector<vm> VMList)
{
	// Placement for offspring
	crossOverCount++;
	Placement child;
	if (VERBOSE)
		cout << "Starting CrossOver on Placements " << par1.PlacementID << " and " << par2.PlacementID << endl;
	unsigned int cutIndex1 = 1 + static_cast<unsigned int>(rand()) % dimPlacement;
	unsigned int cutIndex2;
	if (cutIndex1 > 0)
		cutIndex2 = static_cast<unsigned int>(rand()) % cutIndex1;
	else
		cutIndex2 = 0;
	child.PlacementID = addedPlacementsIteration.size(); // pop.size();
	if (VERBOSE)
	{
		std::cout << "the Placement1 selected for crossOver reaction:";
		par1.printPlacementStructure();
		std::cout << "the Placement2 selected for crossOver reaction:";
		par2.printPlacementStructure();
	}
	// 1-fill the new Placement with X_BS(Placement with best stability level) up to index alphIdx2
	for (unsigned int i = 0; i <= cutIndex2; i++)
	{
		child.structure.push_back(par1.structure[i]);
	}
	// 2-fill the rest of new Placement (>alphIdx2) with the main Placement
	for (unsigned int i = cutIndex2 + 1; i < par2.structure.size(); i++)
	{
		child.structure.push_back(par2.structure[i]);
	}
	if (VERBOSE)
		child.printPlacementStructure();
	// 3-remove repeated indices
	child.putNullValueOnRepeatingVMIndicesHho(cutIndex2); // the repeating indices (from the beginning to apphaIdx2)should be replaced with _ value (null or a specific value eg m where m is the number of VMs and no vm has index m)
	// 4-reinsert missing indices to form a valid structure
	if (VERBOSE)
		child.printPlacementStructure();
	vector<vmPlacement> missingVMs = child.findMissingVMIndicesHho(par1);
	child.insertVMListRandomIntoBubllesHho(missingVMs, par1.structure.size());
	child.calculateFitPlacement(PMList, VMList);
	if (VERBOSE)
	{
		child.printPlacementStructure();
	}
	child.PlacementID = addedPlacementsIteration.size(); // findValidPlacementID();
	updateBestPlacement(child);
	// pop.push_back(New1);
	if (VERBOSE){
		cout << endl
				<< "Placement ID: " << child.PlacementID << endl;
	}
	addedPlacementsIteration.push_back(child);

	/*
	vector<vmPlacement> child_Placement;
	vector<vmPlacement> temp1_Placement;
	vector<vmPlacement> temp2_Placement;
	 int len = par2.structure.size();
	 for(int i = 0;i<len;i++){
		temp1_Placement.push_back(par2.structure[i]);
		temp2_Placement.push_back(structure[i]);
	 }
	 int remaining=len;
	 for(int i = 0;i<len;i++)
	 {
		 // random probability
		 float p = randomNM(0, 100)/100;
		 int selected=rand()%remaining;
		 // if prob is less than 0.45, insert Placement
		 // from parent 1
		 if(p < 0.5){
			int selected_vm=temp1_Placement[selected];
			 child_Placement.push_back(selected_vm);
			 temp1_Placement.erase(temp1_Placement.begin()+selected);
			 temp2_Placement=Placement_remove(temp2_Placement,selected_vm);
			 remaining--;
		 }
		 // if prob is between 0.45 and 0.90, insert
		 // Placement from parent 2
		 else {//if(p < 0.90){
			int selected_vm=temp2_Placement[selected];
			 child_Placement.push_back(selected_vm);
			 temp1_Placement=Placement_remove(temp1_Placement,selected_vm);
			 temp2_Placement.erase(temp2_Placement.begin()+selected);
			 remaining--;
		 }
	 }
	 */
	return child;
};

;

///===========Hho
/*vector<int> createGnome()
{
	int len = num_vm_experiment;
	vector <int> gnome;
	vector<int> list_vm;
	for(int i=0;i<num_vm_experiment;i++) list_vm.push_back(i);
	int num_remaining=num_vm_experiment;
	for(int i=0;i<num_vm_experiment;i++){
		int rand_num=rand()%num_remaining;
		gnome.push_back(list_vm[rand_num]);
		list_vm.erase(list_vm.begin()+rand_num);
		num_remaining--;
		}
	return gnome;
}*/
/*

// Overloading < operator
bool operator<(const Placement &ind1, const Placement &ind2)
{
	return ind1.fit < ind2.fit;
}
float res_wastage_ffd_in_Hho=0;
int num_active_pm_ffd_in_Hho=0;

*/

/*

// Driver code
pair<double,double> main_Hho(vector<pm> pm_set,vector<vm> vm_set)
{
	srand((unsigned)(time(0)));
	for(int i=0;i<num_vm_experiment;i++) Placements[i]=i;
	ofstream Hho_convergence;
	Hho_convergence.open ("conv_Hho.csv");
	Hho_convergence << "iteration,fitness\n";
	// current generation
	int generation = 0;

	vector<Placement> pop;
	bool found = false;
	vector<vm> VMlist2;
	ffd ffd1(pm_set,vm_set,  VMlist2);
	pow_ffd_in_Hho=ffd1.pow_dis;
	// create initial pop
	 cout<<endl<<"gnome:"<<endl;
	for(int i = 0;i<pop_SIZE_Hho;i++)
	{
		vector<int> gnome = create_gnome();
		//for(int j=0;j<num_vm_experiment;j++) cout<<" "<<gnome[j];
	   // cout<<endl;
		pop.push_back(Placement(gnome,pm_set,vm_set));
	}
	int cnt=0;
	int cnt_print=0;
	while(cnt<limit_Hho and FuncEval_Hho<=max_FE_Hho)
	{
		// sort the pop in increasing order of fitness score
		sort(pop.begin(), pop.end());

		// if the Placement having lowest fitness score ie.
		// 0 then we know that we have reached to the target
		// and break the loop

		// Otherwise generate new offsprings for new generation
		vector<Placement> new_generation;

		// Perform Elitism, that mean 10% of fittest pop
		// goes to the next generation
		int s = (10*pop_SIZE_Hho)/100;
		for(int i = 0;i<s;i++)
			new_generation.push_back(pop[i]);

		// From 50% of fittest pop, Placements
		// will mate to produce offspring
		s = (90*pop_SIZE_Hho)/100;
		for(int i = 0;i<s;i++)
		{
		 //   int len = pop.size();
			int r = randomNM(0, 50);
			Placement parent1 = pop[r];
			r = randomNM(0, 50);
			Placement parent2 = pop[r];
			Placement offspring = parent1.mate(parent2, pm_set,  vm_set);
			new_generation.push_back(offspring);
		}
		new_generation=mutate(new_generation);
		pop = new_generation;
	 //   cout<< "Generation: " << generation << "\t";
	//    cout<< "Gnome: "<< endl;
	//    for(int i=0;i<num_vm_experiment;i++) cout<<pop[0].Placement[i] <<" ";
	//    cout<< endl<<"Fitness: "<< pop[0].fitness << "\n";

		generation++;
		cnt++;
		cnt_print++;
		Hho_convergence << generation<<","<<pop[0].fit<<"\n";
		cout<<generation<<"...";
		if(cnt_print==100){
				cout <<endl;
				cnt_print=0;
			   }
	 }
	Hho_convergence.close();
	pair<double,double> ret_pair;
	ret_pair=calculate_stat_Hho( pop[0],pm_set,vm_set);
	 cout<<endl<< "Generation: " << generation << "\t";
   // cout<< "String: "<< pop[0].Placement[0] <<"\t";
	cout<< "Fitness: "<< pop[0].fit << "\n";
	cout<<endl<< "NUM of Functional Evaluations: " << FuncEval_Hho;
	return ret_pair;
}
*/
