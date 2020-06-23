#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include "random.h"
#include "experiment.h"

using namespace std;

Experiment :: Experiment(){}

Experiment :: ~Experiment(){}

void Experiment :: cicleblock(vector<double>& av, Random& rnd, int N, int L){
	vector<double> r;
	double sum;
	int k=0;
	for(int i=0; i<N*L; i++){
		r.push_back(rnd.Rannyu());
	}
	for(int i=0; i<N; i++){
		sum=0;
		for(int j=0; j<L; j++){
			k=j+i*L;
			sum += r[k];
		}
		av.push_back(sum/L);
	}
	return;
}
	
void Experiment :: cicleblock_funz(vector<double>& av, Random& rnd, int N, int L, FunzBase& f){
	vector<double> r;
	double sum;
	int k=0;
	for(int i=0; i<N*L; i++){
		r.push_back(rnd.Rannyu());
	}
	for(int i=0; i<N; i++){
		sum=0;
		for(int j=0; j<L; j++){
			k=j+i*L;
			sum += f.Eval(r[k]);
		}
		av.push_back(sum/L);
	}
	return;
}


