#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "random.h"
#include "Population.h"
#include "Individual.h"

using namespace std;
// constructor
Population :: Population(Random& rnd, int num,vector<double>& x, vector<double>& y){
  Set_Number_ind(num);
  // creando gli individui della popolazione
  for(int i=0; i<Number_ind; i++){
    individuals.push_back(Individual(rnd));
  }
  In_Fitness_order(x,y);
  /*for(int i=0; i<Number_ind; i++){
    cout << individuals[i].Calc_Fitness2_oncirc() << ' ';
    for(int j=0; j<individuals[i].Get_Number_gene(); j++){
      cout << individuals[i].Get_gene(j) << ' ';
    }
    cout << endl;
  }*/
}

// destructor
Population :: ~Population(){}

// settings
void Population :: Set_an_Individual(int a, Individual& Ind){
  for(int i=0; i<individuals[a].Get_Number_gene(); i++){
    individuals[a].Set_Gene(i,Ind.Get_gene(i));
  }
}
// get
double Population :: Average_BHalfpop(vector<double>& x, vector<double>& y){
  double ave=0;
  int half = int((Number_ind/2));
  for(int i=0; i<half; i++){
    ave += individuals[half+i].Calc_Fitness2_oncirc(x,y);
  }
  ave /= half;
  return ave;
}

// methods
int Population :: Selection(Random& rnd){
  double r = pow(rnd.Rannyu(),0.8);
  return int((Number_ind-Number_ind/4)*r);
}

void Population :: New_Generation(vector<Individual>& newgen){
  for(int i=0; i<Number_ind; i++){
    individuals.push_back(newgen[i]);
  }
  return;
}

void Population :: Crossover(Random& rnd, Individual& Son1, Individual& Son2){
  int index1 = Selection(rnd);
  int index2 = Selection(rnd);
  // crossing over point
  point_co = 16;
  // check that index1!=index2
  while(index2==index1){
    index2 = Selection(rnd);
  }
  // Parents
  Individual Parent1 = individuals[index1];
  Individual Parent2 = individuals[index2];
  // Crossover
  int* cross1 = new int[Parent1.Get_Number_gene()-point_co];
  int* cross2 = new int[Parent1.Get_Number_gene()-point_co];
  for(int i=point_co; i<Parent1.Get_Number_gene(); i++){
      cross1[i-point_co] = Parent1.Get_gene(i);
      cross2[i-point_co] = Parent2.Get_gene(i);
  }
  // inserisco in Parent1 in ordine di Parent2 gli elementi di Crossover
  // Parent1--->Son1
  for(int i=0; i<point_co; i++){
    Son1.Set_Gene(i,Parent1.Get_gene(i));
  }
  int ind_gen1 = point_co;
  for(int i=0; i<Parent1.Get_Number_gene(); i++){
    for(int m=0; m<point_co; m++){
      if(Parent2.Get_gene(i)==cross1[m]){
          Son1.Set_Gene(ind_gen1,cross1[m]);
          ind_gen1++;
      }
    }
  }
  // Parent2--->Son2
  for(int i=0; i<point_co; i++){
    Son2.Set_Gene(i,Parent2.Get_gene(i));
  }
  int ind_gen2 = point_co;
  for(int i=0; i<Parent2.Get_Number_gene(); i++){
    for(int m=0; m<point_co; m++){
      if(Parent1.Get_gene(i)==cross2[m]){
          Son2.Set_Gene(ind_gen2,cross2[m]);
          ind_gen2++;
      }
    }
  }
  return;
}

void Population :: Mutation1(Random& rnd, Individual& SonMut){

  // Parent --- > Son
  int ind_parent = int(rnd.Rannyu(0,Number_ind-Number_ind/4));
  Individual ParentMu = individuals[ind_parent];
  for(int i=0; i<ParentMu.Get_Number_gene(); i++){
    SonMut.Set_Gene(i,ParentMu.Get_gene(i));
  }
  // mutation
  // indice di mutazione casuale,
  int index_mu = int(rnd.Rannyu(2,31));
  int app = SonMut.Get_gene(index_mu);
  SonMut.Set_Gene(index_mu,SonMut.Get_gene(index_mu+1));
  SonMut.Set_Gene(index_mu+1,app);
  return;
}

void Population :: Mutation2(Random& rnd, Individual& SonMut){

  // Parent --- > Son
  int ind_parent = int(rnd.Rannyu(0,Number_ind-Number_ind/4));
  Individual ParentMu = individuals[ind_parent];
  // mutation
  SonMut.Set_Gene(0,ParentMu.Get_gene(0));
  for(int i=1; i<(ParentMu.Get_Number_gene()-1); i++){
    SonMut.Set_Gene(i,ParentMu.Get_gene(i+1));
  }
  SonMut.Set_Gene(ParentMu.Get_Number_gene()-1,ParentMu.Get_gene(1));
  return;
}

void Population :: Mutation3(Random& rnd, Individual& SonMut){

  // Parent --- > Son
  int ind_parent = int(rnd.Rannyu(0,Number_ind-Number_ind/4));
  Individual ParentMu = individuals[ind_parent];
  // scelgo un indice tra (1,26) e permuto le prime cinque città da tale indice con le successive in questo modo:
  // 1<-->6; 2<-->7; ...; 5<-->10;
  // mutation or not
  int perm_ind = int(rnd.Rannyu(1,21));
  for(int i=0; i<perm_ind; i++){
    SonMut.Set_Gene(i,ParentMu.Get_gene(i));
  }
  int perm = 0;
  for(int i=perm_ind; i<perm_ind+5; i++){
    SonMut.Set_Gene(i,ParentMu.Get_gene(i+5));
    SonMut.Set_Gene(i+5,ParentMu.Get_gene(i));
    perm++;
  }
  for(int i=perm_ind+10; i<ParentMu.Get_Number_gene(); i++){
    SonMut.Set_Gene(i,ParentMu.Get_gene(i));
  }
  return;
}

void Population :: Mutation4(Random& rnd, Individual& SonMut){

  // Parent --- > Son
  int ind_parent = int(rnd.Rannyu(0,Number_ind-Number_ind/4));
  Individual ParentMu = individuals[ind_parent];
  // scelgo due indici tra (1,31) e ne scambio le città
  // mutation or not
  // scelgo i due indici e mi assicuro che siano diversi
  int index1 = rnd.Rannyu(1,32);
  int index2 = rnd.Rannyu(1,32);
  // check that index1!=index2
  while((index2==index1)){
    index2 = rnd.Rannyu(1,32);
  }
  // Parent --> Son
  for(int i=0; i<ParentMu.Get_Number_gene(); i++){
    SonMut.Set_Gene(i,ParentMu.Get_gene(i));
  }
  // mutation
  int appoggio = SonMut.Get_gene(index1);
  SonMut.Set_Gene(index1,SonMut.Get_gene(index2));
  SonMut.Set_Gene(index2,appoggio);
  return;
}

void Population :: In_Fitness_order(vector<double>& x, vector<double>& y){
  vector<Individual> F_order;
  double fittest=0, last_fittest=0;
  int index_f=0, order=0;
  for(int i=0; i<Number_ind; i++){
    if(individuals[i].Calc_Fitness2_oncirc(x,y)>fittest){
      fittest = individuals[i].Calc_Fitness2_oncirc(x,y);
      index_f = i;
    }
  }
  F_order.push_back(individuals[index_f]);
  order++;
  while(order<Number_ind){
    last_fittest=fittest;
    fittest = 0;
    for(int i=0; i<Number_ind; i++){
      double fit = individuals[i].Calc_Fitness2_oncirc(x,y);
      if(fit > fittest && fit < last_fittest){
        fittest = fit;
        index_f = i;
      }
    }
    F_order.push_back(individuals[index_f]);
    order++;
  }
  for(int i=0; i<Number_ind; i++){
    individuals[i] = F_order[i];
  }
  F_order.clear();
  /*for(int i=0; i<Number_ind; i++){
    cout << individuals[i].Calc_Fitness2_oncirc(x,y) << ' ';
    for(int j=0; j<individuals[i].Get_Number_gene(); j++){
      cout << individuals[i].Get_gene(j) << ' ';
    }
    cout << endl;
  }*/
  return;
}
