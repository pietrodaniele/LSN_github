#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "random.h"
#include "Population.h"
#include "Individual.h"
#include "funzioni.h"
#include <mpi.h>

using namespace std;
// quantità utili
int n_ind[2] = {2000,5000}, N_gnrts = 0, early_stop=0, f=0, pie=0, Ng=0;
int generations[2] = {1200,800};
double prop_mut = 0.09, prop_cross = 0.7, last_betsfitness = 0;
vector<double> x, y;
// operatori
Random rnd , rnd2;
// funzioni

int main (int argc, char *argv[]){
  if(argc!=2){
    cerr << "./main.exe 0 or 1" << endl;
    cerr << "0 ---> 32 cities randomly placed on a circumference" << endl;
    cerr << "1 ---> 32 cities randomly placed inside a square" << endl;
    exit(1);
  }
  int config = atoi(argv[1]);
  int indiv = n_ind[config];
  // setto il generatore rnd uguale per tutti i continenti in modo da avere tutti la stessa config di città
  createrandom(rnd2);
  createrandom(rnd);
  Cities_config(rnd,config,x,y);
  // Setto le variabili per il calcolo parallelo
  int size, rank, tag1=1, tag2=2;
  MPI_Status status1, status2;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int* send_array = new int[32];
  int* rec_array = new int[32];
  string ranks[4] = {"rank0","rank1","rank2","rank3"};
  int sr_0=1, sr_1=0, sr_2=3, sr_3=2; // defisco gli scambi
  // setto il generatore del random, tutte le grandezze utili
  createrandom_mpi(rnd,rank);
  // apro i file dove andrò a scrivere il miglior cammino e la fitness sulla metà della bestpop
  ofstream Write_Best, Write_Ave;
  if(config==0){
    Write_Best.open("data/Best_indiv_"+ranks[rank]+"_circ.0");
    Write_Ave.open("data/Ave_half_"+ranks[rank]+"_circ.0");
  }else{
    Write_Best.open("data/Best_indiv_"+ranks[rank]+"_square.0");
    Write_Ave.open("data/Ave_half_"+ranks[rank]+"_square.0");
  }
	// creo la generazione 0, di partenza
  cout << "Every genenation has " << indiv << " individuals." << endl;
  Population* population = new Population(rnd,indiv,x,y);
  population->In_Fitness_order(x,y);
  cout << "Genenation " << N_gnrts <<", rank = " << rank << " ---> Fittest = " << population->Get_BestFitnees(x,y) << endl;
  Write_Ave << N_gnrts << " " << population->Average_BHalfpop(x,y) << endl;
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
    Write_Ave << N_gnrts << " " << population->Average_BHalfpop(x,y) << endl;
    cout << "Genenation " << N_gnrts <<", rank = " << rank << " ---> Fittest = " << population->Get_BestFitnees(x,y) << endl;
    if(Ng==50){
      // estraggo 4 numeri casuali tutti diversi ad esempio (2,0,3,1), così trovo quali continenti si scambiano i best path
      // in questo caso avrei 2<-->0 3<-->1.
      /*int rnd_rank1, rnd_rank2, rnd_rank3, rnd_rank4;
      rnd_exc(rnd2,rnd_rank1,rnd_rank2,rnd_rank3,rnd_rank4);
      int exc[4] = {rnd_rank1, rnd_rank2, rnd_rank3, rnd_rank4};*/
      if (rank == 0){
        // carico l'array da inviare
        for(int i=0; i<32; i++){
            send_array[i] = population->Get_an_Individual(indiv-1).Get_gene(i);
        }
        MPI_Send(send_array,32,MPI_INT,sr_0,tag1,MPI_COMM_WORLD);
        MPI_Recv(rec_array,32,MPI_INT,sr_0,tag2,MPI_COMM_WORLD,&status1);
        population->Set_an_Individual_array(indiv-1,rec_array);
        population->In_Fitness_order(x,y);
      }
      if (rank == 1){
        for(int i=0; i<32; i++){
            send_array[i] = population->Get_an_Individual(indiv-1).Get_gene(i);
        }
        MPI_Recv(rec_array,32,MPI_INT,sr_1,tag1,MPI_COMM_WORLD,&status2);
        MPI_Send(send_array,32,MPI_INT,sr_1,tag2,MPI_COMM_WORLD);
        population->Set_an_Individual_array(indiv-1,rec_array);
        population->In_Fitness_order(x,y);
      }
      if (rank == 2){
        // carico l'array da inviare
        for(int i=0; i<32; i++){
            send_array[i] = population->Get_an_Individual(indiv-1).Get_gene(i);
        }
        MPI_Send(send_array,32,MPI_INT,sr_2,tag1,MPI_COMM_WORLD);
        MPI_Recv(rec_array,32,MPI_INT,sr_2,tag2,MPI_COMM_WORLD,&status1);
        population->Set_an_Individual_array(indiv-1,rec_array);
        population->In_Fitness_order(x,y);
      }
      if (rank == 3){
        for(int i=0; i<32; i++){
            send_array[i] = population->Get_an_Individual(indiv-1).Get_gene(i);
        }
        MPI_Recv(rec_array,32,MPI_INT,sr_3,tag1,MPI_COMM_WORLD,&status2);
        MPI_Send(send_array,32,MPI_INT,sr_3,tag2,MPI_COMM_WORLD);
        population->Set_an_Individual_array(indiv-1,rec_array);
        population->In_Fitness_order(x,y);
      }
      cout << " Array sent" << endl;
      Ng=0;
    }
    pie++;
    Ng++;
    if(pie==(int)generations[config]/2){
      sr_0=3;
      sr_1=2;
      sr_2=1;
      sr_3=0;
    }
  }
  for(int i=0; i<32; i++){
    int city = population->Get_an_Individual(indiv-1).Get_gene(i);
    Write_Best << city << " " << x[city-1] << " " << y[city-1] << endl;
  }
  Write_Best.close();
  Write_Ave.close();
  MPI_Finalize();
return 0;
}
