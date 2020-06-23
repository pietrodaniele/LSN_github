#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "random.h"
#include "metropolis.h"
#include "FunzBase3D.h"
#include "hydro_100.h"
#include "hydro_210.h"
#include "funzioni.h"

using namespace std;
// quantità utili
double a_0 = 0.0529;
double x_0=a_0, y_0=a_0, z_0=a_0, delta100=1.17*a_0, delta210=2.8*a_0, delta100_gauss=0.72*a_0, delta210_gauss=2*a_0;
int n = 1000000, L=10000;
vector<double> x, y, z, r, r_2, sum_prog, sum2_prog, error;
// operatori
Metropolis metro;
Random rnd;
// funzioni
Hydro_100 hydro100;
Hydro_210 hydro210;

int main (){
	// setto il generatore del random, il mio cube lattice
    createrandom(rnd);
	// setto ora tutti i vari parametri, tipo delta e punto di partenza
	metro.SetDelta(delta100);
	metro.SetStart(x_0,y_0,z_0);
	// ora faccio partire l'algoritmo di metropolis che mi stimerà le distribuzioni di probabilità
	// n=1 l=0 m=0
	cout << "WaveFunction n=1 l=0 m=0 (1s)" << endl;
	cout << "Uniform distribution" << endl;
	metro.Algoritm_unif(n,hydro100,rnd,x,y,z);
	metro.coord_r(r,r_2,x,y,z,L);
	metro.accumulation(sum_prog,r);
	metro.accumulation(sum2_prog,r_2);
	metro.errorprog(sum_prog,sum2_prog,error);
	write_100(x,y,z,sum_prog,error);
	clean(x,y,z,r,r_2,sum_prog,sum2_prog,error);
	// studio per il caso 100 anche una distribuzione gaussiana
	cout << "Gaussian distribution" << endl;
	metro.SetDelta(delta100_gauss);
	metro.Algoritm_gauss(n,hydro100,rnd,x,y,z);
	metro.coord_r(r,r_2,x,y,z,L);
	metro.accumulation(sum_prog,r);
	metro.accumulation(sum2_prog,r_2);
	metro.errorprog(sum_prog,sum2_prog,error);
	write_100gauss(x,y,z,sum_prog,error);
	clean(x,y,z,r,r_2,sum_prog,sum2_prog,error);
	// risetto i parametri
	cout << "WaveFunction n=2 l=1 m=0 (2s)" << endl;
	cout << "Uniform distribution" << endl;
	metro.SetDelta(delta210);
	metro.SetStart(5*x_0,5*y_0,5*z_0);
	metro.Algoritm_unif(n,hydro210,rnd,x,y,z);
	metro.coord_r(r,r_2,x,y,z,L);
	metro.accumulation(sum_prog,r);
	metro.accumulation(sum2_prog,r_2);
	metro.errorprog(sum_prog,sum2_prog,error);
	write_210(x,y,z,sum_prog,error);
	clean(x,y,z,r,r_2,sum_prog,sum2_prog,error);
  // studio per il caso 210 anche una distribuzione gaussiana
  cout << "Gaussian distribution" << endl;
  metro.SetDelta(delta210_gauss);
  metro.Algoritm_gauss(n,hydro210,rnd,x,y,z);
  metro.coord_r(r,r_2,x,y,z,L);
  metro.accumulation(sum_prog,r);
  metro.accumulation(sum2_prog,r_2);
  metro.errorprog(sum_prog,sum2_prog,error);
  write_210gauss(x,y,z,sum_prog,error);
  clean(x,y,z,r,r_2,sum_prog,sum2_prog,error);
return 0;
}
