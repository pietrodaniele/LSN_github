#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include "integranda.h"
#include "experiment.h"
#include "funzioni.h"

using namespace std;

Random rnd;
Experiment lab;
Integranda funz;
vector<double> av, av2, err_prog, av_S, err_S_prog, sum_prog, sum_S_prog, sum2_prog;
int M=10000, N=100, L=M/N, a=0, b=1;

int main (int argc, char *argv[]){
	// preparo il generatore random
	createrandom(rnd);
	/* ora procedo a calcolare l'integrale con uniform e sampling, ovviamente il sampling serve solo per una migliore stima dell'incertezza.
	Calcolando con lo stesso punto random riesco ad avere un confronto a parità di dati delle incertezze. Per la distribuzione non unif è 
	stato usato usato un semplice andamento lineare p(x)∝(1-x) --> p(x)=2(1-x) e pertanto g'(x)=g(x)/(2*(1-x)). Così ad ogni stima dell'in=
	grale ho le due incertezze.*/
	// calcolo l'integrale per g(x) e g'(x)
	lab.cicleblock_integral(av,av_S,rnd,N,L,a,b,funz);
	// ora calcolo val_quad per unif accumulo e trovo l'errore. Poi svuoto i vector che voglio riutilizzare
	lab.valoriquad(av,av2);
	lab.accumulation(sum_prog,av);
	lab.accumulation(sum2_prog,av2);
	lab.errorprog(sum_prog,sum2_prog,err_prog);
	clean(av2,sum2_prog);
	// ora calcolo val_quad per non unif accumulo e trovo l'errore.
	lab.valoriquad(av_S,av2);
	lab.accumulation(sum_S_prog,av_S);
	lab.accumulation(sum2_prog,av2);
	lab.errorprog(sum_S_prog,sum2_prog,err_S_prog);
	// ora scrivo il file dove ho la stima dell'integrale e le due incertezze
	write(sum_prog,err_prog,err_S_prog);
	return 0;
}
