#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "random.h"
#include "ising_1D.h"

using namespace std;

void createrandom(Random&);

double minimum(double, double);

int Pbc(int, int);  //Algorithm for periodic boundary conditions

void calculation(Ising_1D&, ofstream&, double*,  double*, double*, double*);

double error(double*, double*,int);
