#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "FunzBase3D.h"
#include "funzioni.h"
#include "random.h"
#include "metropolis.h"

using namespace std;

Metropolis :: Metropolis(){}

Metropolis :: ~Metropolis(){}



void Metropolis ::  SetDelta(double a){
	delta=a;
	return;
}

void Metropolis ::  SetStart(double a, double b, double c){
	x_0=a;
	y_0=b;
	z_0=c;
	return;
}

void Metropolis :: Algoritm_unif(int n, FunzBase3D& f, Random& rnd, vector<double>& x, vector<double>& y, vector<double>& z){
	double alfa, x_try, y_try, z_try, acc=0;
	// setto l'origine
	x.push_back(x_0);
	y.push_back(y_0);
	z.push_back(z_0);
	for(int i=1; i<n; i++){
		// fornisce il punto try
		x_try = x[i-1]+rnd.Rannyu(-delta,delta);
		y_try = y[i-1]+rnd.Rannyu(-delta,delta);
		z_try = z[i-1]+rnd.Rannyu(-delta,delta);
		alfa = minimum(1,f.Eval(x_try,y_try,z_try)/f.Eval(x[i-1],y[i-1],z[i-1]));
		accettanza_uniform(alfa,x,y,z,rnd,acc,i,x_try,y_try,z_try);
	}
	// stampo ora una stima dell'accettanza
	cout << "L'accettanza è del " << acc*100/n << "%" << endl;
	return;
}

void Metropolis :: Algoritm_gauss(int n, FunzBase3D& f, Random& rnd, vector<double>& x, vector<double>& y, vector<double>& z){
	double alfa, x_try, y_try, z_try, acc=0;
	// setto l'origine
	x.push_back(x_0);
	y.push_back(y_0);
	z.push_back(z_0);
	for(int i=1; i<n; i++){
		// fornisce il punto try
		x_try = x[i-1]+rnd.Gauss(0,delta);
		y_try = y[i-1]+rnd.Gauss(0,delta);
		z_try = z[i-1]+rnd.Gauss(0,delta);
		alfa = minimum(1,f.Eval(x_try,y_try,z_try)/f.Eval(x[i-1],y[i-1],z[i-1]));
		accettanza_uniform(alfa,x,y,z,rnd,acc,i,x_try,y_try,z_try);
	}
	// stampo ora una stima dell'accettanza
	cout << "L'accettanza è del " << acc*100/n << "%" << endl;
	return;
}

void Metropolis :: coord_r(vector<double>& r, vector<double>& r_2, vector<double>& x, vector<double>& y, vector<double>& z, int L){
	double sum;
	unsigned int n = x.size()/L;
	for(unsigned int i=0; i<n; i++){
		sum=0;
		for(unsigned int j=i*L; j<i*L+L; j++){
			sum+=pow(x[j]*x[j]+y[j]*y[j]+z[j]*z[j],0.5);
			// cout  << i <<" --- " << j << endl;
		}
		r.push_back(sum/L);
		r_2.push_back(r[i]*r[i]);
		if(i%10==0 && i!=0){
			cout <<"Number of block = "<< i << endl;
		}
	}
	cout << "Number of block = "<< r.size() << endl << endl;
	return;
}

void Metropolis :: accumulation(vector<double>& sum_prog, vector<double>& av){
	for(unsigned int i=0; i<av.size(); i++){
		sum_prog.push_back(0);
	}
	// cout << "azzerato" << endl;
	for(unsigned int i=0; i<av.size(); i++){
		for(unsigned int j=0; j<i+1; j++){
			sum_prog[i]+=av[j];
		}
		sum_prog[i]=sum_prog[i]/(i+1);
		if(i%100==0){
		}
	}
	// cout << "Accumulation done" << endl << endl;
	return;
}

void Metropolis :: errorprog(vector<double>& sum_prog, vector<double>& sum2_prog, vector<double>& error_prog){
	for(unsigned int i=0; i<sum_prog.size(); i++){
		if(i==0){
			error_prog.push_back(0);
		 }else{
			error_prog.push_back(sqrt((sum2_prog[i]-(sum_prog[i]*sum_prog[i]))/i));
		}
	}
	return;
}
