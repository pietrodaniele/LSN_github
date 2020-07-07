#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "random.h"
#include "SA.h"

using namespace std;
// constructor
SimAnneling :: SimAnneling(){}
// destructor
SimAnneling :: ~SimAnneling(){}
// settings

// methods
void SimAnneling :: Anneal(Individual& Ind, Random& rnd, vector<double>& x, vector<double>& y, int m){
  // writing metric during the SA
  ofstream Write_L2;
  string type[2] = {"circ","square"};
  Write_L2.open("data/L2_"+type[m]+".dat");
  // set insta Temp and the evolution counter
  Temp = start_Temp;
  int ev_count=1;
  do{
    Individual New_Ind(rnd,32);
    New_Ind.Set_Genis(Ind.Get_Genis());
    New_Ind.Evolve(rnd);
    double p = exp(-(New_Ind.Calc_Fitness2(x,y)-Ind.Calc_Fitness2(x,y))/Temp); // K_b is equal 1
    if(rnd.Rannyu() < p){
      Ind.Set_Genis(New_Ind.Get_Genis());
    }
    // writing L2 at the {ev_count} evolution
    Write_L2 << Temp << " " << ev_count << " " << Ind.Calc_Fitness2(x,y) << endl;
    Temp = Temp - deltaTemp;
    ev_count++;
  }while(Temp > absTemp+deltaTemp);
  Write_L2.close();
}
// città distribuite equi-spaziate su una circo unitaria
