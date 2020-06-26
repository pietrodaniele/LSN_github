#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "FunzBase.h"
#include "random.h"


using namespace std;

void createrandom(Random&);

void accettanza_uniform(double, vector<double>&, Random&, double&, int, double);

void clean(vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&);

double minimum(double, double);

void write_try(vector<double>&, vector<double>&);

void norm(FunzBase&, Random& );
