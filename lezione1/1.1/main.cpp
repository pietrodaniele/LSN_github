#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include "random.h"
#include "FunzBase.h"
#include "FunzSigma.h"
#include "experiment.h"
#include "funzioni.h"


using namespace std;

Random rnd;
Experiment lab;
vector<double> av, av2, sum_prog, sum2_prog, error_prog, chiquadro;
int M=1000000, N=100, L=M/N, Q=200, n=10000;
FunzSigma sigma;
 
int main (){
	createrandom(rnd);
	// prima parte
	lab.cicleblock(av,av2,rnd,N,L); // primo ciclo sul numero dei blocchi
	lab.accumulation(sum_prog,av,N); // accumullo per av
	lab.accumulation(sum2_prog,av2,N); // accumulo per av2
	lab.errorprog(sum_prog,sum2_prog,error_prog,N); // calcolo le varianze dalla differenza del valore quadratrico medio e del valore medio al quadrato
	write(sum_prog,error_prog); // scrivo x sum_prog e error_prog su un file da carica sui notebook
	// pulisco i miei vector per poterli riutilizzarli
	clean(av,av2,sum_prog,sum2_prog,error_prog);
	// seconda parte
	lab.cicleblock_funz(av,av2,rnd,N,L,sigma); // primo ciclo sul numero dei blocchi però ora con una funzione
	lab.accumulation(sum_prog,av,N); // accumullo per av
	lab.accumulation(sum2_prog,av2,N); // accumulo per av2
	lab.errorprog(sum_prog,sum2_prog,error_prog,N); // calcolo le varianze dalla differenza del valore quadratrico medio e del valore medio al quadrato
	write_sigma(sum_prog,error_prog); // scrivo x sum_prog e error_prog su un file da carica sui notebook
	// terza parte ---> calcolo il chiquadro
	lab.chiquadro(rnd,Q,n,chiquadro);
	write_chiquadro(chiquadro);
	return 0;
}


