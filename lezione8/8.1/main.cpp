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
vector<double> x, x2, sum_prog, sum2_prog, error;
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
	// ora faccio partire l'algoritmo di metropolis che mi stimerà le distribuzioni di probabilità
	// mu = 0.5 and sigma = 1
	cout << "WaveFunction_try" << endl;
	cout << "Uniform distribution" << endl;
  cout << "Starting point x = " << s << endl;
  cout << "Delta = " << D << endl;
  Try.Set_Mu(0.74599);
  Try.Set_Sigma(0.518303);
  cout << "Setting mu and sigma:" << endl;
  cout << "Mu = " << Try.Get_Mu() << endl;
  cout << "Sigma = " << Try.Get_Sigma() << endl;
	metro.Algoritm_unif(N*L,Try,rnd,x);
  metro.Algoritm_unif_Try(N*L,Try,rnd,x2);
  write_try(x,x2);
  cout << "Metropolis done..." << endl;
  metro.cicleblock_integral(sum_prog,x,N,L,Try);
  norm(Try,rnd);
	//clean(x,x2,sum_prog,sum2_prog,error)
	// studio per il caso 100 anche una distribuzione gaussiana
return 0;
}
