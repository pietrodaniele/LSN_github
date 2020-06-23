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
int indiv = 1000, N_gnrts = 0, early_stop=0, f=0, pie=0;
double prop_mut = 0.01, prop_cross = 0.9, last_betsfitness = 0;
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
	// setto il generatore del random, tutte le grandezze utili
  createrandom(rnd);
  Cities_config(rnd,config,x,y);
  // apro i file dove andrò a scrivere il miglior cammino e la fitness sulla metà della bestpop
  ofstream Write_Best, Write_Ave;
  if(config==0){
    Write_Best.open("data/Best_indiv_circ.0");
    Write_Ave.open("data/Ave_half_circ.0");
  }else{
    Write_Best.open("data/Best_indiv_square.0");
    Write_Ave.open("data/Ave_half_square.0");
  }
	// creo la generazione 0, di partenza
  cout << "Every genenation has " << indiv << " individuals." << endl;
  Population* population = new Population(rnd,indiv,x,y);
  cout << "Genenation " << N_gnrts << " ---> Fittest = " << population->Get_BestFitnees(x,y) << endl;
  Write_Ave << N_gnrts << " " << population->Average_BHalfpop(x,y) << endl;
  Population* newgen = new Population(rnd,indiv,x,y);
  while(pie<1000){
    Population* appog = new Population(rnd,indiv,x,y);
    if(f==0){
      appog = population;
      f++;
    }else{
      appog = newgen;
    }
    int newgen_count = 0;
    while(newgen_count<indiv-1){
      Individual Son1(rnd), Son2(rnd), SonMut(rnd);
      /*if(newgen_count<indiv){
        cout
        appog->Mutation1(rnd,SonMut,prop_mut);
        newgen->Set_an_Individual(newgen_count,SonMut);
        newgen_count++;
      }else if(newgen_count<indiv){
        appog->Mutation2(rnd,SonMut,prop_mut);
        newgen->Set_an_Individual(newgen_count,SonMut);
        newgen_count++;
      }else if(newgen_count<indiv){
        appog->Mutation3(rnd,SonMut,prop_mut);
        newgen->Set_an_Individual(newgen_count,SonMut);
        newgen_count++;
      }else if(newgen_count<indiv){
        appog->Mutation4(rnd,SonMut,prop_mut);
        newgen->Set_an_Individual(newgen_count,SonMut);
        newgen_count++;
      }else if(newgen_count<indiv){
        cout << "Cross" << endl;
        appog->Crossover(rnd,Son1,Son2, prop_cross);
        if(newgen_count==indiv-1){
          newgen->Set_an_Individual(newgen_count,Son1);
          newgen_count++;
        }else{
          newgen->Set_an_Individual(newgen_count,Son1);
          newgen_count++;
          newgen->Set_an_Individual(newgen_count,Son2);
          newgen_count++;
        }
      }
    }
    N_gnrts++;
    //cout << N_gnrts << endl;
    //Print(population);
    newgen->In_Fitness_order(x,y);
    cout << "Genenation " << N_gnrts << " ---> Fittest = " << newgen->Get_BestFitnees(x,y) << endl;
    Write_Ave << N_gnrts << " " << newgen->Average_BHalfpop(x,y) << endl;
    double deltafit = abs(newgen->Get_BestFitnees(x,y)-last_betsfitness);
    if(deltafit <= 0.001){
      early_stop++;
    }else{
      early_stop=0;
    }
    last_betsfitness = newgen->Get_BestFitnees(x,y);
    pie++;
  }
  for(int i=0; i<32; i++){
    int city = newgen->Get_an_Individual(indiv-1).Get_gene(i);
    Write_Best << city << " " << x[city-1] << " " << y[city-1] << endl;
  }
  Write_Best.close();
  Write_Ave.close();
return 0;
}
