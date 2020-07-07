#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "random.h"
#include "Population.h"
#include "Individual.h"

using namespace std;

void createrandom(Random&);

void createrandom_mpi(Random&, int);

void Cities_config(Random&, int, vector<double>&, vector<double>&);

//void New_Generation(Population*, Population*);

void Print(Population*,vector<double>&, vector<double>&);

Individual Mutation(Population*, Random&, vector<double>&, vector<double>&, double);

bool Is_New_Son(Individual, vector<Individual>&);

void rnd_exc(Random&, int&, int&, int&, int&);
//void Write_BHalf_pop_fit(Population*);

//void Write_Bestfitness(Population*,vector<double>& x, vector<double>& y);
