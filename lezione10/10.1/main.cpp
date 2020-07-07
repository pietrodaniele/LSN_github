#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "random.h"
#include "SA.h"
#include "Individual.h"
#include "funzioni.h"
#include "mpi.h"

using namespace std;
// quantità utili
double S_Temp, d_Temp, absTemp;
int Ng=32;
vector<double> x, y;
// operatori
Random rnd;
SimAnneling SAnn;
// funzioni

int main (int argc, char const *argv[]){
  if(argc!=2){
    cerr << "./main.exe 0 or 1" << endl;
    cerr << "0 ---> 32 cities randomly placed on a circumference" << endl;
    cerr << "1 ---> 32 cities randomly placed inside a square" << endl;
    exit(1);
  }
  if(atoi(argv[1])==0){
    S_Temp=10;
    d_Temp=0.0001;
    absTemp=0.0001;
  }else{
    S_Temp=1;
    d_Temp=0.00001;
    absTemp=0.00001;
  }
  // setting the rnd gen
  createrandom(rnd);
  // creating cities config
  Cities_config(rnd,atoi(argv[1]),Ng,x,y);
  // setting SA parameters
  SAnn.Set_start_Temp(S_Temp);
  SAnn.Set_dT(d_Temp);
  SAnn.Set_absTemp(absTemp);
  // creating Individual
  Individual Ind(rnd,Ng);
  // show parameters
  print_parameters(S_Temp,d_Temp,absTemp,Ng,atoi(argv[1]));
  // SA
  SAnn.Anneal(Ind,rnd,x,y,atoi(argv[1])); // mi fornisce il migliore percorso
  // Write the best path
  Ind.Write_path(atoi(argv[1]),x,y);
  cout << "Best path at Temp = " << absTemp << " :" << endl;
  Ind.Print();
  cout << "L^2 = " << Ind.Calc_Fitness2(x,y) << endl;
return 0;
}
