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
};

void print_parameters(unsigned int M, unsigned int n, unsigned int l, unsigned int t_steps,\
                      double S0, double T, double K, double r, double sigma){
  cout << "MC steps M = "<< M << endl;
  cout << "# of blocks N = " << n << endl;
  cout << "# of MC steps in each block = " << l << endl;
  cout << endl;
  cout << "---------------------------------------" << endl;
  cout << endl;
  cout << "Market's parameters:" << endl;
  cout << "		- asset price at t=0: S(0)=" << S0 << endl;
  cout << "		- delivery time: T=" << T << endl;
  cout << "		- strike price: K="<< K << endl;
  cout << "		- risk-free interest rate: r=" << r << endl;
  cout << "		- volatility: sigma=" << sigma << endl;
  cout << "                - time steps for discret calc=" << t_steps << endl;
  cout << endl;
return;
}
