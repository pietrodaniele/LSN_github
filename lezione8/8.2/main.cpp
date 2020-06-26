#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "random.h"
#include "metropolis.h"
#include "FunzBase.h"
#include "funz_try.h"
#include "funzioni.h"

using namespace std;
// quantità utili
int N = 100, L=10000;
vector<double> x, x2, av;
// operatori
Metropolis metro;
Random rnd;
// funzioni
Funz_Try Try;

int main (){
	// setto il generatore del random, il mio cube lattice
  createrandom(rnd);
	// setto ora tutti i vari parametri, tipo delta e punto di partenza
  double s=3, D=2;
	metro.SetDelta(D);
	metro.SetStart(s);
  // trovo i migliori parametri all'interno di mu=[0.8.0.85], sigma=[0.6,0.65]
  double best_mu=0, best_sigma=0, best_H=0;
  vector<double> H_t;
  cout << "Searching best params" << endl;
  for(double a=0.8; a<=0.85; a+=0.002){ // ciclo su mu
    for(double b=0.6; b<=0.65; b+=0.002){ // ciclo su sigma
      Try.Set_Mu(a);
      Try.Set_Sigma(b);
      // studio le distribuzioni
      metro.Algoritm_unif(N*L,Try,rnd,x);
      metro.Algoritm_unif_Try(N*L,Try,rnd,x2);
      // metroplis <H>_t
      metro.cicleblock_integral(H_t,x,N,L,Try);
      if(H_t[N-1]<=best_H){
        best_mu=a;
        best_sigma=b;
        best_H=H_t[N-1];
        H_t.clear();
      }else{
        H_t.clear();
      }
      x.clear();
      x2.clear();
    }
  }
  cout << "Best params: Mu = " << best_mu << "; Sigma = " << best_sigma << endl;
	// ora faccio partire l'algoritmo di metropolis che mi stimerà le distribuzioni di probabilità
	// mu = 0.5 and sigma = 1
	cout << "WaveFunction_try" << endl;
	cout << "Uniform distribution" << endl;
  cout << "Starting point x = " << s << endl;
  cout << "Delta = " << D << endl;
  Try.Set_Mu(best_mu);
  Try.Set_Sigma(best_sigma);
  // studio le distribuzioni
  metro.Algoritm_unif(N*L,Try,rnd,x);
  metro.Algoritm_unif_Try(N*L,Try,rnd,x2);
  //Try.Set_Mu(0.811);
  //Try.Set_Sigma(0.621);
  cout << "Setting mu and sigma:" << endl;
  cout << "Mu = " << Try.Get_Mu() << endl;
  cout << "Sigma = " << Try.Get_Sigma() << endl;
  write_try(x,x2);
  //ene(Try);
  cout << "Metropolis done..." << endl;
  metro.cicleblock_integral(av,x,N,L,Try);
  cout << "<H>_t value = " << av[N-1] << endl;
  metro.calc_err(av);
  // scrivo gli andamenti del phi e phi^2
  norm(Try,rnd);
	//clean(x,x2,sum_prog,sum2_prog,error)
	// studio per il caso 100 anche una distribuzione gaussiana
return 0;
}
