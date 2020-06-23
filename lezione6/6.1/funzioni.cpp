#include "funzioni.h"
#include <fstream>
#include <string>

using namespace std;

void createrandom(Random& rnd){
   int seed[6];
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

double minimum(double a, double b){
	if( a <= b )
		return a;
	else
		return b;
};

int Pbc(int i, int nspin){
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
};

void calculation(Ising_1D& ising, ofstream& Write, double* ene, double* heat, double* magn, double* susc){
	// definisco array necessari
	int N = ising.Get_nblock();
	double* av2 = new double[N];
	double* sum = new double[N];
	double* sum2 = new double[N];
	double* err = new double[N];
	// calcolo la somma prog e l'errore
	ising.av2(ene,av2);
	ising.accumulation(ene,sum,av2,sum2,err);
	Write << sum[N-1] << " " << err[N-1] << " ";
	ising.av2(heat,av2);
	ising.accumulation(heat,sum,av2,sum2,err);
	Write << sum[N-1] << " " << err[N-1] << " ";
	ising.av2(magn,av2);
	ising.accumulation(magn,sum,av2,sum2,err);
	Write << sum[N-1] << " " << err[N-1]  << " ";
	ising.av2(susc,av2);
	ising.accumulation(susc,sum,av2,sum2,err);
	Write << sum[N-1] << " " << err[N-1]  << endl;
	return;
};


double error(double* av, double* av2, int n){
    if(n==0)
        return 0;
    else
        return sqrt((av2[n] - av[n]*av[n])/n);
};
