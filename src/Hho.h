/*
 * Hho.h
 *
 *  Created on: Jan 15, 2024
 *      Author: mahdi
 */

#ifndef HHO_H_
#define HHO_H_



#include <iostream>
#include <vector>
#include <cmath>
#include <stdio.h>	/* printf, scanf, puts, NULL */
#include <stdlib.h> /* srand, rand */
#include <time.h>	/* time */
#include <iterator>
#include <vector>
#include <algorithm>
#include "defs.h"
#include "mapping.h"

#include <random>
#include <iomanip>
#include <string>
#include <map>
#include <bits/stdc++.h>
#include <fstream>

#include "ffd.h"
#include "Placement.h"

using namespace std;

// Valid Genes
// int Placements[num_vm_experiment];



class Hho
{
private:
public:
	Hho();
	virtual ~Hho();
	Hho(int epoch, int pop_size, int num_dims);
	vector<Placement> pop;
	//void Hhoin(int pop_size);
	//vector<Placement> population;
    //Placement bestSolution;
    default_random_engine generator;
    int energy;
	unsigned int epoch;
	unsigned int currentIteration;
	unsigned int pop_size;
	bool sort_flag;
	unsigned int dimPlacement;
	float bestFit;
	float worstFit;
	unsigned int numFuncEval_Hho;
	unsigned int crossOverCount, mutationCount;
	const double EPSILON = 1e-8;
	std::vector<Placement> addedPlacementsIteration;
	std::vector<Placement> toBeRemovedPlacementsIteration;
	Placement bestPlacement;
	Placement worstPlacement;
	Placement bestEnrichmentBound; // X_SB
	Placement PlacementBestStabLevel;
	int mappingSol[num_vm_experiment];
	int numActive_pm;
	mapping Solution;
	float bestPow_dis;
	float bestResourceWastage;
	float bestStabLevel;
	float enrichmentBound;
	void run( vector<vm> vms,  vector<pm> servers);
	void scoutPhase( vector<vm> vms,  vector<pm> servers);
	void explorePhase( vector<vm> vms,  vector<pm> servers);
	void exploitPhase( vector<vm> vms,  vector<pm> servers);
	void attackPhase( vector<vm> vms,  vector<pm> servers);
	void evolve(vector<pm> PMList, vector<vm> VMList, std::ofstream &fileName, unsigned int iterationForAvg);
	void initialize(vector<vm> VMList);
	void calculateFitnessPop(vector<pm> PMList, vector<vm> VMList); // calculate Placements' fitness and find best and worst Placements and save them in bestPlacement and worstPlacement Placements
	void sortPopFitnessDecreasing();								// calculate Placements' fitness and find best and worst Placements and save them in bestPlacement and worstPlacement Placements
	vector<Placement> sortPopByFitnessDecreasingly();
	unsigned int findPlacementIdInPop(Placement thePlacement);
	void printPop();
	void printPlacementVector(vector<Placement> vectPlacements);
	void printFitPop();
	unsigned int findValidPlacementID(); // search into pop and find free Placement IDs and return a valid one
	void updateBestPlacement(Placement theNewPlacement);
	void updateWorstPlacement(Placement theNewPlacement);
	void modifyPlacementIndices(); // to correct Placements IDs and give them IDs from 0 up
	void setBestAndWorstPlacements();
	void positronReactionOld(Placement thePlacement, vector<pm> PMList, vector<vm> VMList);
	vector<Placement> mutate(vector<Placement> popPlacements);
	Placement explore(Placement par1, Placement par2, vector<pm> PMList, vector<vm> VMList);
};




#endif /* Hho_H_ */
