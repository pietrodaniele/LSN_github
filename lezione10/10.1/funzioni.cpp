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

void Cities_config(Random& rnd, int config, int Ng,vector<double>& x, vector<double>& y){
  ofstream Write;
  if(config==0){
    Write.open("data/Cities_config_circ.0");
    for(int i=0; i<Ng; i++){
      double teta = rnd.Rannyu(0,2*M_PI);
      x.push_back(cos(teta));
      y.push_back(sin(teta));
      Write << i+1 << " " << x[i] << " " << y[i] << endl;
    }
  }else if(config==1){
    Write.open("data/Cities_config_square.0");
    for(int i=0; i<Ng; i++){
      x.push_back(rnd.Rannyu());
      y.push_back(rnd.Rannyu());
      Write << i+1 << " " << x[i] << " " << y[i] << endl;
    }
  }
  Write.close();
};

void print_parameters(double st, double dt, double at, int Nc, int m){
  string type[2] = {"Circ","Square"};
  cout << "System parameters:" << endl;
  cout << "   - Start Temp = " << st << endl;
  cout << "   - Cooling rate dt = " << dt << endl;
  cout << "   - Absolute Temp = " << at << endl;
  cout << "   - Number of cities = " << Nc << endl;
  cout << "   - Cities config = " << type[m] << endl;
  cout << "   - Metric = L^2 " << endl;
  cout << endl;
};
