#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "random.h"
#include "FunzBase.h"
#include "FunzExp.h"
#include "FunzLnz.h"
#include "experiment.h"
#include "funzioni.h"

using namespace std;

Random rnd;
Experiment lab;
FunzExp esp;
FunzLnz lnz;
vector<double> v1, v2, v10, v100;
int N=10000;
 
int main (int argc, char *argv[]){
//////////////////////PREPARO IL GEN RANDOM///////////////////////
    createrandom(rnd);
	// UNIFORME
    lab.cicleblock(v1,rnd,N,1);
    lab.cicleblock(v2,rnd,N,2);
    lab.cicleblock(v10,rnd,N,10);
    lab.cicleblock(v100,rnd,N,100);
    write_unif(v1,v2,v10,v100);
    clean(v1,v2,v10,v100);
    // EXP
    lab.cicleblock_funz(v1,rnd,N,1,esp);
    lab.cicleblock_funz(v2,rnd,N,2,esp);
    lab.cicleblock_funz(v10,rnd,N,10,esp);
    lab.cicleblock_funz(v100,rnd,N,100,esp);
    write_exp(v1,v2,v10,v100);
    clean(v1,v2,v10,v100);    
	// CAUCHY LORENTZ
	lab.cicleblock_funz(v1,rnd,N,1,lnz);
	lab.cicleblock_funz(v2,rnd,N,2,lnz);
    lab.cicleblock_funz(v10,rnd,N,10,lnz);
    lab.cicleblock_funz(v100,rnd,N,100,lnz);
    write_lnz(v1,v2,v10,v100);
    clean(v1,v2,v10,v100);
	return 0;
}

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//implemento ora la distribuzione exp e cauchy-lorentz
/*	for(int i=0;i<P;i++){
		exp_x[i]=(-1/lambda)*log(1-rnd.Rannyu()); // esponenziale
		cl_x[i]=T*tan(M_PI*(rnd.Rannyu()-1/2));   // cauchy-lorentz
	}
	ofstream fout("dati.dat");
		for(int i=0; i<P ; i++){
			fout<<exp_x[i]<<" "<<cl_x[i]<<endl;
	}
	fout.close();*/
