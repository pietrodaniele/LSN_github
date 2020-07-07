#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "random.h"
#include "Population.h"
#include "Individual.h"
#include "funzioni.h"

using namespace std;
// quantità utili
int n_ind[2] = {2000,5000}, N_gnrts = 0, early_stop=0, f=0, pie=0;
int generations[2] = {1200,800};
double prop_mut = 0.09, prop_cross = 0.7, last_betsfitness = 0;
vector<double> x, y;
// operatori
Random rnd;
// funzioni

int main (int argc, char const *argv[]){
  if(argc!=2){
    cerr << "./main.exe 0 or 1" << endl;
    cerr << "0 ---> 32 cities randomly placed on a circumference" << endl;
    cerr << "1 ---> 32 cities randomly placed inside a square" << endl;
    exit(1);
  }
  int config = atoi(argv[1]);
  int indiv = n_ind[config];
	// setto il generatore del random, tutte le grandezze utili
  createrandom(rnd);
  Cities_config(rnd,config,x,y);
  // apro i file dove andrò a scrivere il miglior cammino e la fitness sulla metà della bestpop
  ofstream Write_Best, Write_Ave, Write_L2;
  if(config==0){
    Write_Best.open("data/Best_indiv_circ.0");
    Write_Ave.open("data/Ave_half_circ.0");
    Write_L2.open("data/L2_circ.0");
  }else{
    Write_Best.open("data/Best_indiv_square.0");
    Write_Ave.open("data/Ave_half_square.0");
    Write_L2.open("data/L2_square.0");
  }
	// creo la generazione 0, di partenza
  cout << "Every genenation has " << indiv << " individuals." << endl;
  Population* population = new Population(rnd,indiv,x,y);
  population->In_Fitness_order(x,y);
  cout << "Genenation " << N_gnrts << " ---> Fittest = " << population->Get_BestFitnees(x,y) << endl;
  Write_Ave << N_gnrts << " " << population->Average_BHalfpop(x,y) << endl;
  Write_L2 << N_gnrts << " " <<  population->Get_BestFitnees(x,y) << endl;
  vector<Individual> newgen;
  while(pie<generations[config]){
    int newgen_count=0, m=0;
    //riempio la popolazione con i primi N migliori individui non mutati
    while(newgen_count<indiv/4){
      Individual Ind = population->Get_an_Individual(indiv-1-m);
      m++;
      bool Truth = Is_New_Son(Ind,newgen);
      while(Truth!=true){
        Ind = population->Get_an_Individual(indiv-1-m);
        m++;
        Truth = Is_New_Son(Ind,newgen);
      }
      newgen.push_back(Ind);
      newgen_count++;
    }
    for(int i=indiv/4; i<indiv; i++){
      Individual Prova = Mutation(population,rnd,x,y,prop_mut);
      bool Truth = Is_New_Son(Prova,newgen);
      while(Truth!=true){
        Prova = Mutation(population,rnd,x,y,prop_mut);
        Truth = Is_New_Son(Prova,newgen);
      }
        newgen.push_back(Prova);
    }
    population->Clean_population();
    population->New_Generation(newgen);
    newgen.clear();
    population->In_Fitness_order(x,y);
    N_gnrts++;
    Write_L2 << N_gnrts << " " <<  population->Get_BestFitnees(x,y) << endl;
    Write_Ave << N_gnrts << " " << population->Average_BHalfpop(x,y) << endl;
    cout << "Genenation " << N_gnrts << " ---> Fittest = " << population->Get_BestFitnees(x,y) << endl;
    pie++;
  }
  for(int i=0; i<32; i++){
    int city = population->Get_an_Individual(indiv-1).Get_gene(i);
    Write_Best << city << " " << x[city-1] << " " << y[city-1] << endl;
  }
  Write_Best.close();
  Write_Ave.close();
  Write_L2.close();
return 0;
}
