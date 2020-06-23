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
	
void Experiment :: cicleblock_integral(vector<double>& av, vector<double>& avS, Random& rnd, int N, int L, int a, int b, FunzBase& f){
	double sum, sumS, x;
	for(int i=0; i<N; i++){
		sum=0;
		sumS=0;
		for(int j=0; j<L; j++){
			x=rnd.Rannyu(a,b);
			sum += f.Eval(x);
			sumS +=f.Eval_sampling(x);
		}
		av.push_back(sum/L);
		avS.push_back(sumS/L);
	}
	return;
}

void Experiment :: valoriquad(vector<double>& av, vector<double>& av2){
	for(unsigned int i=0; i<av.size(); i++){
		av2.push_back(av[i]*av[i]);
	}
	return;
}


void Experiment :: accumulation(vector<double>& sum_prog, vector<double>& av){
	for(unsigned int i=0; i<av.size(); i++){
		sum_prog.push_back(0);
	}
	for(unsigned int i=0; i<av.size(); i++){
		for(unsigned int j=0; j<i+1; j++){
			sum_prog[i]+=av[j];
		}
		sum_prog[i]=sum_prog[i]/(i+1);
	}
	return;	
}

void Experiment :: errorprog(vector<double>& sum_prog, vector<double>& sum2_prog, vector<double>& error_prog){
	for(unsigned int i=0; i<sum_prog.size(); i++){
		if(i==0){
			error_prog.push_back(0);
		 }else{
			error_prog.push_back(sqrt((sum2_prog[i]-(sum_prog[i]*sum_prog[i]))/i));
		}
	}
	return;
}
