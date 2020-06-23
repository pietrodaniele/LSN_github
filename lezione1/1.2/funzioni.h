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

void write_unif(vector<double> v1, vector<double> v2, vector<double> v10, vector<double> v100){
	ofstream WriteData;
	WriteData.open("data_unif.dat");
	for(unsigned int i=0; i<v1.size(); i++){
		WriteData<<v1[i]<<" "<<v2[i]<<" "<<v10[i]<<" "<<v100[i]<<endl;
	}
	WriteData.close();
	return;
};

void write_exp(vector<double> v1, vector<double> v2, vector<double> v10, vector<double> v100){
	ofstream WriteData;
	WriteData.open("data_exp.dat");
	for(unsigned int i=0; i<v1.size(); i++){
		WriteData<<v1[i]<<" "<<v2[i]<<" "<<v10[i]<<" "<<v100[i]<<endl;
	}
	WriteData.close();
	return;
};

void write_lnz(vector<double> v1, vector<double> v2, vector<double> v10, vector<double> v100){
	ofstream WriteData;
	WriteData.open("data_lnz.dat");
	for(unsigned int i=0; i<v1.size(); i++){
		WriteData<<v1[i]<<" "<<v2[i]<<" "<<v10[i]<<" "<<v100[i]<<endl;
	}
	WriteData.close();
	return;
};

void clean(vector<double>& av, vector<double>& av2, vector<double>& sum_prog, vector<double>& sum2_prog){
	av.clear();
	av2.clear();
	sum_prog.clear();
	sum2_prog.clear();
	return;	
};

#endif // __funzioni__
