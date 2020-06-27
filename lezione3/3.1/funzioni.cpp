#include "funzioni.h"
#include <cmath>

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

void write(vector<double>& sum_prog, vector<double>& error_prog){
	ofstream WriteData;
	WriteData.open("data.dat");
	for(unsigned int i=0; i<sum_prog.size(); i++){
		WriteData<<i<<" "<<sum_prog[i]<<" "<<error_prog[i]<<endl;
	}
	WriteData.close();
	return;
};

double d(double S, double t, double r, double sigma, double K, double T) {
	return (log(S/K)+r+sigma*sigma*(T-t)/2.)/(sigma*sqrt(T-t));
};

double N(double x) {
	return 0.5*(1+erf(x/sqrt(2)));
}
