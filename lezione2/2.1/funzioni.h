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

void clean(vector<double>& av2, vector<double>& sum2_prog){
	av2.clear();
	sum2_prog.clear();
	return;	
};

void write(vector<double> v1, vector<double> v2, vector<double> v3){
	ofstream WriteData;
	WriteData.open("data.dat");
	for(unsigned int i=0; i<v1.size(); i++){
		WriteData<<i<<" "<<v1[i]<<" "<<v2[i]<<" "<<v3[i]<<endl;
	}
	WriteData.close();
	return;
};

#endif // __funzioni__
