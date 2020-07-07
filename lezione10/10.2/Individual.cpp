#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "random.h"
#include "Individual.h"

using namespace std;
// constructor
Individual :: Individual(Random& rnd){
  Genis.push_back(1);
  int p;
  for(int i=1; i<Number_gene; i++){
    int test=0;
    while(test!=1){
      p = int(rnd.Rannyu(2,33));
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
}
// destructor
Individual :: ~Individual(){}
// settings

// methods
// città distribuite equi-spaziate su una circo unitaria

double Individual :: Calc_Fitness2_oncirc(vector<double>& x, vector<double>& y){
    double Fit2 = 0;
    for(unsigned int i=0; i<(Genis.size()-1); i++){
      double dist2 = pow(x[Genis[i+1]-1]-x[Genis[i]-1],2) + pow(y[Genis[i+1]-1]-y[Genis[i]-1],2);
      Fit2 += dist2;
    }
    Fit2 += pow(x[Genis[Genis.size()-1]-1]-x[Genis[0]-1],2) + pow(y[Genis[Genis.size()-1]-1]-y[Genis[0]-1],2);
    return Fit2;
}

double Individual :: Calc_Fitness2_oncirc_print(vector<double>& x, vector<double>& y){
    double Fit2 = 0;
    for(unsigned int i=0; i<(Genis.size()-1); i++){
      cout << Genis[i] << endl;
      double dist2 = pow(x[Genis[i+1]-1]-x[Genis[i]-1],2) + pow(y[Genis[i+1]-1]-y[Genis[i]-1],2);
      Fit2 += dist2;
      cout << Fit2 << endl;
    }
    Fit2 += pow(x[Genis[Genis.size()-1]-1]-x[Genis[0]-1],2) + pow(y[Genis[Genis.size()-1]-1]-y[Genis[0]-1],2);
    return Fit2;
}
