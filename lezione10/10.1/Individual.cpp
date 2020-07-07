#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "random.h"
#include "Individual.h"

using namespace std;
// constructor
Individual :: Individual(Random& rnd, int Ng){
  Number_gene = Ng;
  Genis.push_back(1);
  int p;
  for(int i=1; i<Number_gene; i++){
    int test=0;
    while(test!=1){
      p = int(rnd.Rannyu(2,Number_gene+1));
      for(unsigned int i=0; i<Genis.size(); i++){
        if(p==Genis[i]){
          test=2;
        }
      }
      if(test==2){
        test=0;
      }else{
        test=1;
      }
    }
    Genis.push_back(p);
  }
  //cout << "Starting path:" << endl;
  //Print();
}
// destructor
Individual :: ~Individual(){}
// settings
void Individual :: Set_Genis(vector<int> a){
  for(int i=0; i<Number_gene; i++){
    Genis[i]=a[i];
  }
}
// methods
double Individual :: Calc_Fitness2(vector<double>& x, vector<double>& y){
    double Fit2 = 0;
    for(unsigned int i=0; i<(Genis.size()-1); i++){
      double dist2 = pow(x[Genis[i+1]-1]-x[Genis[i]-1],2) + pow(y[Genis[i+1]-1]-y[Genis[i]-1],2);
      Fit2 += dist2;
    }
    Fit2 += pow(x[Genis[Genis.size()-1]-1]-x[Genis[0]-1],2) + pow(y[Genis[Genis.size()-1]-1]-y[Genis[0]-1],2);
    return Fit2;
}
//
void Individual :: Evolve(Random& rnd){
  double r = rnd.Rannyu();
  if(r < 0.25){
    int i = (int)rnd.Rannyu(1,32), j = (int)rnd.Rannyu(1,32);
    Mutation1(i,j);
  }else if(r>=0.25 && r<0.5){
    int start = (int)rnd.Rannyu(2,21), n_cities = (int)rnd.Rannyu(1,12);
    Mutation2(n_cities,start);
  }else if(r>=0.5 && r<0.75){
    int i = (int)rnd.Rannyu(2,21), nblk = (int)rnd.Rannyu(1,12);
    Mutation3(i,nblk);
  }else if(r>=0.75){
    int k = (int)rnd.Rannyu(2,21), n_inv = (int)rnd.Rannyu(1,12);
    Mutation4(k,n_inv);
  }
};
// scambio
void Individual :: Mutation1(int i, int j){
  int app = Genis[i];
  Genis[i] = Genis[j];
  Genis[j] = app;
  return;
};
// shift 1 of {n_cities} cities starting from index start
void Individual :: Mutation2(int n_cities, int start){
  if(start+n_cities>=Number_gene){
    cerr << "Error Mut2" << endl;
    exit(1);
  }
  int app[n_cities],  m=0, g=0;
  for(int i=start; i<start+n_cities; i++){
    app[m] = Genis[i];
    m++;
  }
  for(int j=start; j<start+n_cities-1; j++){
    Genis[j+1] = app[g];
    g++;
  }
  Genis[start]=app[g];
}
// permutation (i,i+n_block)
void Individual :: Mutation3(int i, int N_block){
  int t, n=N_block/2;
  for(int j=i; j<(int)(i+N_block)/2; j++){
    t = Genis[j];
    Genis[j] = Genis[j+n];
    Genis[j+n] = t;
  }
}
// invert
void Individual :: Mutation4(int index, int n_inv){
  int t;
	for (int i=0; i<(int)(n_inv/2); i++) {
	   t = Genis[index+i];
	   Genis[i+index] = Genis[index+n_inv-1-i];
	   Genis[index+n_inv-1-i] = t;
		}
}
// write Genis
void Individual :: Write_path(int o, vector<double>& x, vector<double>& y){
  ofstream Write;
  string type[2] = {"circ","square"};
  Write.open("data/best_path_"+type[o]+".dat");
  for(int i=0; i<Number_gene; i++){
    int city = Genis[i];
    Write << city << " " << x[city-1] << " " << y[city-1] << endl;
  }
}
// print
void Individual :: Print(){
  for(int i=0; i<Number_gene; i++){
    cout << Genis[i] << " ";
  }
  cout << endl;
}
