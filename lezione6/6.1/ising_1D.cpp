#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "funzioni.h"
#include "random.h"
#include "ising_1D.h"

using namespace std;

Ising_1D :: Ising_1D(){}

Ising_1D :: ~Ising_1D(){}

void Ising_1D :: Input(Random& rnd){
	ifstream ReadInput;

  	cout << "Classic 1D Ising model             " << endl;
  	cout << "Monte Carlo simulation             " << endl << endl;
  	cout << "Nearest neighbour interaction      " << endl << endl;
  	cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  	cout << "The program uses k_B=1 and mu_B=1 units " << endl;
  	// Read input informations
  	// Temperatures
  	ReadInput.open("input.dat");
  	ReadInput >> temp_low;
  	beta_low = 1.0/temp_low;
  	cout << "Temperature low = " << temp_low << endl;
  	ReadInput >> temp_high;
  	beta_high = 1.0/temp_high;
  	cout << "Temperature high = " << temp_high << endl;
  	// Spin e system's values
  	ReadInput >> nspin;
  	cout << "Number of spins = " << nspin << endl;
  	ReadInput >> J;
  	cout << "Exchange interaction = " << J << endl;
  	ReadInput >> h;
  	cout << "External field = " << h << endl << endl;
  	// Use Gibbs or Metropolis
  	ReadInput >> metro; // if=1 Metropolis else Gibbs
  	// Setting blocks and MC steps
  	ReadInput >> nblock;
  	ReadInput >> nstep;
	// Measure between temp_low e temp_high
	ReadInput >> camp;
	cout << "I campionamenti tra " << temp_low << " e " << temp_high << " sono: " << camp << endl << endl;
	// metro or gibbs
  	if(metro==1) cout << "The program perform Metropolis moves" << endl;
  	else cout << "The program perform Gibbs moves" << endl;
  	cout << "Number of blocks = " << nblock << endl;
  	cout << "Number of steps in one block = " << nstep << endl << endl;
  	ReadInput.close();
  	// Initial configuration
  	for (int i=0; i<nspin; ++i){
  		if(rnd.Rannyu() >= 0.5) s[i] = 1;
  	    else s[i] = -1;
  	}
  	/*/Evaluate energy etc. of the initial configuration
  	  Measure();

  	//Print initial values for the potential energy and virial
  	  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;*/
};

void Ising_1D :: Move(double beta, Random& rnd){
	int o;
	double energy_old, energy_new, p, A;
	double energy_up, energy_down;
	// Faccio nspin tentativi di girare lo spin. Con metro non ho l'accettanza al 100%, con gibbs si.
	for(int i=0; i<nspin; ++i){
  		//Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    	o = (int)(rnd.Rannyu()*nspin);
    	//Metropolis
    	if(metro==1){
			energy_old = -J * s[o] * ( s[Pbc(o-1,nspin)] + s[Pbc(o+1,nspin)] ) - h * s[o];
			energy_new = -J * (-s[o]) * ( s[Pbc(o-1,nspin)] + s[Pbc(o+1,nspin)] ) - h *(-s[o]);
			p = exp(-beta*(energy_new-energy_old));
			A = minimum(1.0,p);
			double r = rnd.Rannyu();
			if(r < A)
				s[o] = -s[o];
    	}
    	//Gibbs sampling
    	else{
				energy_up = -J * s[o] * ( s[Pbc(o-1,nspin)] + s[Pbc(o+1,nspin)] ) - h * s[o];
				energy_down = -J * (-s[o]) * ( s[Pbc(o-1,nspin)] + s[Pbc(o+1,nspin)] ) - h * (-s[o]);
				p = 1/(1.+exp(-beta*(energy_up-energy_down)));
				if(p > rnd.Rannyu()){
					s[o] = -s[o];
				}
    	}
  	}
};

void Ising_1D :: Measure(){
	//double e_var[nspin], var=0, e_mean=0;
	double ene_2 = 0;
	// muovo dopo la misura in quanto ho già la config iniziale
	for(int j=0; j<nspin; j++){
		double ene = -J * s[j] * s[Pbc(j+1,nspin)] - 0.5 * h * (s[j] + s[Pbc(j+1,nspin)]);
		u += ene;
		/*e_var[j]=ene;
		e_mean += ene;*/
		ene_2 += pow(ene,2);
		magn += s[j];
		double Sp = 0.0, Uj = 0.0;
		for(int b=0; b<nspin; b++){
			Sp += s[b];
		}
		chi += s[j]*Sp;
	}
	u2 += ene_2/(nspin*nstep);
	/*e_mean /= nspin;
	for(int j=0; j<nspin; j++){
		var += (e_var[j]-e_mean)*(e_var[j]-e_mean);
	}
	u2 += var/nspin;*/
};

void Ising_1D :: av2(double* av, double* av2){
	for(int i=0; i<nblock; i++){
		av2[i]=av[i]*av[i];
	}
};

void Ising_1D :: accumulation(double* av, double* sum_prog, double* av2, double* sum2_prog, double* err_prog){
	for(int i=0; i<nblock; i++){
		sum_prog[i]=0;
		sum2_prog[i]=0; // azzero i vettori sum prog per sicurezza
	    for(int j=0; j<i+1; j++){
	        sum_prog[i] += av[j];
	        sum2_prog[i] += av2[j];
		}
		sum_prog[i]/=(i+1);
    	sum2_prog[i]/=(i+1);
    	err_prog[i] = error(sum_prog,sum2_prog,i);
	}
};

void Ising_1D :: Auto_eq(Random& rnd, double temp, double sigma){
	// ciclo di equilibrazione
	cout << endl;
	cout << "Equilibrating the configuration at the temperature " << temp << " ..." << endl << endl;
	double Av = 0;
	unsigned int MCsteps = 1, early_stop=0;
	int k=1;
	do{
		u=0;
		Measure();
		double old_Av = Av;
		Av = Av*(nspin)*(k-1) + u;
		Av /= (nspin*k);
		if((old_Av-Av)<sigma){
			early_stop++;
		}else{
			early_stop=0;
		}
		MCsteps++;
		k++;
		Move(1/temp,rnd);
	}while(early_stop!=10 && MCsteps<=100000);
	cout << "The config. at the temperature " << temp << " is equilibrated with a tolerance="<<sigma<<" after " << MCsteps << " steps." << endl<< endl;

	return;
};
