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

void Experiment :: cicleblock(vector<double>& av, vector<double>& av2, Random rnd, int N, int L){
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
		av2.push_back(av[i]*av[i]);
	}
	return;
}

void Experiment :: accumulation(vector<double>& sum_prog, vector<double>& av, int N){
	for(int i=0; i<N; i++){
		sum_prog.push_back(0);
	}
	for(int i=0; i<N; i++){
		for(int j=0; j<i+1; j++){
			sum_prog[i]+=av[j];
		}
		sum_prog[i]=sum_prog[i]/(i+1);
	}
	return;	
}

void Experiment :: errorprog(vector<double>& sum_prog, vector<double>& sum2_prog, vector<double>& error_prog, int N){
	for(int i=0; i<N; i++){
		if(i==0){
			error_prog.push_back(0);
		 }else{
			error_prog.push_back(sqrt((sum2_prog[i]-(sum_prog[i]*sum_prog[i]))/i));
		}
	}
	return;
}
	



