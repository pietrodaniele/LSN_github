#include "funzioni.h"
#include "random.h"
#include <string>

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

void accettanza_uniform(double alfa, vector<double>& x, Random& rnd, double& acc, int i, double x_try){
	double r=rnd.Rannyu();
	if(r<=alfa){
		x.push_back(x_try);
		acc++;
	}else{
		x.push_back(x[i-1]);
	}
	return;
};

void clean(vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& v4, vector<double>& v5){
	x.clear();
	y.clear();
	z.clear();
	v4.clear();
	v5.clear();
	return;
};

double minimum(double a, double b){
	if(a<=b){
		return a;
	}else{
		return b;
	}
};

void write_try(vector<double>& x, vector<double>& x2){
	ofstream write;
	write.open("data/funz_try.dat");
	for(unsigned int h=0; h<x.size(); h++){
		write << x[h] << " " << x2[h] << endl;
	}
  write.close();
  return;
}

void norm(FunzBase& f, Random& rnd){
  ofstream Norm, Norm2;
  Norm.open("data/Funz.norm");
  Norm2.open("data/Funz2.norm");
  for(double i=-7.5 ;i<=7.5; i+=0.05){
    Norm << f.Eval(i) << endl;
    Norm2 << f.Eval2(i) << endl;

  }
  Norm.close();
  Norm2.close();
  return;
}
