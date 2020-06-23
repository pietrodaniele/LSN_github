#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "FunzBase.h"
#include "funzioni.h"
#include "random.h"
#include "metropolis.h"

using namespace std;
// constructor
Metropolis :: Metropolis(){}
// destructor
Metropolis :: ~Metropolis(){}
// settings
void Metropolis ::  SetDelta(double a){
	delta=a;
	return;
}

void Metropolis ::  SetStart(double a){
	x_0=a;
	return;
}
// methods
void Metropolis :: Algoritm_unif(int n, FunzBase& f, Random& rnd, vector<double>& x){
	double alfa, x_try, acc=0;
	// setto l'origine
	x.push_back(x_0);
	for(int i=1; i<n; i++){
		// fornisce il punto try
		x_try = x[i-1]+rnd.Rannyu(-delta,delta);
		alfa = minimum(1,f.Eval2(x_try)/f.Eval2(x[i-1]));
		accettanza_uniform(alfa,x,rnd,acc,i,x_try);
	}
	// stampo ora una stima dell'accettanza
	cout << "L'accettanza è del " << acc*100/n << "%" << endl;
	return;
}

void Metropolis :: Algoritm_unif_Try(int n, FunzBase& f, Random& rnd, vector<double>& x){
	double alfa, x_try, acc=0;
	// setto l'origine
	x.push_back(x_0);
	for(int i=1; i<n; i++){
		// fornisce il punto try
		x_try = x[i-1]+rnd.Rannyu(-delta,delta);
		alfa = minimum(1,f.Eval(x_try)/f.Eval(x[i-1]));
		accettanza_uniform(alfa,x,rnd,acc,i,x_try);
	}
	// stampo ora una stima dell'accettanza
	cout << "L'accettanza è del " << acc*100/n << "%" << endl;
	return;
}

void Metropolis :: cicleblock_integral(vector<double>& av, vector<double>& x, int N, int L,FunzBase& f){
	double sum, sum_norm=0;
	for(unsigned int i=0;i<x.size();i++){
		sum_norm += f.Eval2(x[i]);
	}
	double norm = sum_norm/x.size();
	//cout << norm << endl;
	for(int i=0; i<N; i++){
		sum=0;
		for(int j=0; j<L; j++){
			int k = i*L+j;
			double H = -0.5*f.EvalD2(x[k])+(pow(x[k],4)-(5./2.)*pow(x[k],2))*f.Eval(x[k]);
			//double H = (pow(x[k],4)-(5./2.)*pow(x[k],2))*f.Eval(x[k]);
			//sum += f.Eval2(x[k])*H/(f.Eval(x[k]));
			sum += H/(f.Eval(x[k]));

		}
		av.push_back(sum/(L));
		//cout << av[i] << endl;
	}
	return;
}

void Metropolis :: calc_err(vector<double>& av){
	double* av2 = new double[av.size()];
	double* sum = new double[av.size()];
	double* sum2 = new double[av.size()];
	double* err = new double[av.size()];
	ofstream Write;
	Write.open("data/output.blk");
	// init sums array
	for(unsigned int i=0; i<av.size(); i++){
		sum[i]=0;
		sum2[i]=0;
	}
	// av2
	for(unsigned int i=0; i<av.size(); i++){
		av2[i]=av[i]*av[i];
	}
	// acc
	for(unsigned int i=0; i<av.size(); i++){
		for(unsigned int j=0; j<i+1; j++){
			sum[i]+=av[j];
			sum2[i]+=av2[j];
		}
		sum[i]/=(i+1);
    sum2[i]/=(i+1);
		if(i==0){
			err[i]=0;
		}else{
			err[i] = sqrt((sum2[i] - sum[i]*sum[i])/i);
		}
		Write << sum[i] << " " << err[i] << endl;
	}
	Write.close();
}
