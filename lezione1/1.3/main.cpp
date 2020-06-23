#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "random.h"
#include "FunzBase.h"
#include "experiment.h"
#include "funzioni.h"

using namespace std;

Random rnd;
Experiment lab;
vector<double> av, av2, sum_prog, sum2_prog, error_prog;
int M=1000000, N=100, L=M/N;
double l=0.5, d=2;
 
int main (int argc, char *argv[]){
	// setto il generatore del random
    	createrandom(rnd);
    	// ora faccio i cicli sui blocchi e accumulo
   	lab.cicleblock_buffon(av,av2,rnd,N,L,l,d); // nel metodo ho imposto che l'intervallo diviso è [0,10)
   	lab.accumulation(sum_prog,av);
    	lab.accumulation(sum2_prog,av2);
	// ora faccio gli errori
	lab.errorprog(sum_prog,sum2_prog,error_prog);
	// infine scrivo il file
	write(sum_prog,error_prog);
return 0;
}
