#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "random.h"

using namespace std;

void createrandom(Random&);

void accettanza_uniform(double, vector<double>&, vector<double>&, vector<double>&, Random&, double&, int, double, double, double);

void write_100(vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&);

void write_100gauss(vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&);

void write_210(vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&);

void write_210gauss(vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&);

void clean(vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&);

double minimum(double, double);
