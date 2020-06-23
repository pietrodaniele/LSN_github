#ifndef __funzioni__
#define __funzioni__

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
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

/*double error(vector<double>& AV, vector<double>& AV2, int n){
	if(n==0)
		return 0;
	else
		return sqrt((AV2[n]-(AV[n]*AV[n]))/n);
};*/

void write(vector<double> sum_prog, vector<double> error_prog){
	ofstream WriteData;
	WriteData.open("data.dat");
	for(unsigned int i=0; i<sum_prog.size(); i++){
		WriteData<<i<<" "<<sum_prog[i]<<" "<<error_prog[i]<<endl;
	}
	WriteData.close();
	return;
};


#endif // __funzioni__
