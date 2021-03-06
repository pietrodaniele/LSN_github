#ifndef __funzioni__
#define __funzioni__

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "random.h"

using namespace std;

void createrandom(Random&);

void write(vector<double>&, vector<double>&);

double d(double, double, double, double, double, double);

double N(double);

void print_parameters(unsigned int, unsigned int, unsigned int, unsigned int,\
                      double, double, double, double, double);

#endif // __funzioni__
