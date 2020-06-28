#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "experiment.h"

using namespace std;

Experiment :: Experiment(){}

Experiment :: ~Experiment(){}

void Experiment :: cicleblock(double *av, double *av2, string *type, int nmis, int L, int k, string s){
	double sum, x;
	ifstream Read;
	Read.open("data/output_"+type[k]+"_"+s+".dat");
	for(int i=0; i<nmis/L; i++){
		sum=0;
		for(int j=0; j<L; j++){
			Read >> x;
			sum += x;
		}
		av[i]=(sum/L);
		av2[i]=(av[i]*av[i]);;
	}
	Read.close();
	return;
}

void Experiment :: accumulation(double *sum_prog, double *av, int nmis){
	for(int i=0; i<nmis; i++){
		sum_prog[i]=0;
	}
	for(int i=0; i<nmis; i++){
		for(int j=0; j<i+1; j++){
			sum_prog[i]+=av[j];
		}
		sum_prog[i]=sum_prog[i]/(i+1);
	}
	return;
}

void Experiment :: errorprog(double* sum_prog, double *sum2_prog, double *error_prog, int nmis){
	for(int i=0; i<nmis; i++){
		if(i==0){
			error_prog[i]=0;
		 }else{
			error_prog[i]=(sqrt((sum2_prog[i]-(sum_prog[i]*sum_prog[i]))/i));
		}
	}
	return;
}
