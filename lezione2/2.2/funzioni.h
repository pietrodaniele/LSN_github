#ifndef __funzioni__
#define __funzioni__

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;

void createrandom(Random& rnd){
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;
   rnd.SaveSeed();
};

void write_lattice(vector<double> v1, vector<double> v2){
	ofstream WriteData;
	WriteData.open("data_lattice.dat");
	for(unsigned int i=0; i<v1.size(); i++){
		WriteData<<v1[i]<<" "<<v2[i]<<endl;
	}
	WriteData.close();
	return;
};

void write_continuum(vector<double> v1, vector<double> v2){
	ofstream WriteData;
	WriteData.open("data_continuum.dat");
	for(unsigned int i=0; i<v1.size(); i++){
		WriteData<<v1[i]<<" "<<v2[i]<<endl;
	}
	WriteData.close();
	return;
};

void clean3(vector<double>& av, vector<double>& av2, vector<double>& sum_prog){
	av.clear();
	av2.clear();
	sum_prog.clear();
	return;	
};

void clean2(vector<double>& av, vector<double>& av2){
	av.clear();
	av2.clear();
	return;	
};

void error(vector<double>& r, vector<double>& error){
	double r_2[r.size()], sum[r.size()], sum2[r.size()];
	for(unsigned int i=0; i<r.size(); i++){
			r_2[i]=(r[i]*r[i]);
			sum[i]=0;
			sum2[i]=0;
	}
	for(unsigned int i=0; i<r.size(); i++){
		for(unsigned int j=0; j<i+1; j++){
			sum[i]+=r[j];
			sum2[i]+=r_2[j];
		}
		sum[i] = sum[i]/(i+1);
		sum2[i] = sum2[i]/(i+1);
		if(i==0){
			error.push_back(0);
		}else{
			error.push_back(sqrt((sum2[i]-(sum[i]*sum[i]))/i));
		}
	}
	return;
};

void r2_mean_lattice(int N, int M, Random& rnd, RandomWalk& drunk, vector<double>& r, vector<double>& x, vector<double>& y, vector<double>& z){
	for(int i=0; i<M; i++){
		drunk.Walk_cube(rnd,N,x,y,z);
	    if(i==0){
	    	for(int j=0; j<N; j++){
	    	r.push_back(x[j]*x[j]+y[j]*y[j]+z[j]*z[j]);
	    	}
	    }else{
	    	for(int j=0; j<N; j++){
	    		r[j] = (r[j]*i + (x[j]*x[j]+y[j]*y[j]+z[j]*z[j]))/(i+1);
	    	}
	    }
	    clean3(x,y,z);
	}
	return;
};

void r2_mean_continuum(int N, int M, Random& rnd, RandomWalk& drunk, vector<double>& r, vector<double>& x, vector<double>& y, vector<double>& z){
	for(int i=0; i<M; i++){
		drunk.Walk_continuum(rnd,N,x,y,z);
	    if(i==0){
	    	for(int j=0; j<N; j++){
	    	r.push_back(x[j]*x[j]+y[j]*y[j]+z[j]*z[j]);
	    	}
	    }else{
	    	for(int j=0; j<N; j++){
	    		r[j] = (r[j]*i + (x[j]*x[j]+y[j]*y[j]+z[j]*z[j]))/(i+1);
	    	}
	    }
	    clean3(x,y,z);
	}
	return;
};

#endif // __funzioni__
