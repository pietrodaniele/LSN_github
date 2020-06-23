#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "random.h"
#include "ising_1D.h"
#include "funzioni.h"

using namespace std;
// parameters
Random rnd;
Ising_1D ising;


int main(){
	createrandom(rnd);
	ising.Input(rnd);
	// file dove saranno salvati i valori delle grandezze e le loro incertezze
	ofstream Write;
	if(ising.Get_metro()==1){
		if(ising.Get_h()==0)
			Write.open("values_metro.0");
		else
			Write.open("values_metro_h.0");
	}else{
		if(ising.Get_h()==0)
			Write.open("values_gibbs.0");
		else
			Write.open("values_gibbs_h.0");}
	// definisco tutti i vettori che mi servono ene, ene2, heat,
	double* ene = new double[ising.Get_nblock()];
	double* ene2 = new double[ising.Get_nblock()];
	double* heat = new double[ising.Get_nblock()];
	double* magnet = new double[ising.Get_nblock()];
	double* susc = new double[ising.Get_nblock()];
	// intervallo di temp tra un campionamento e l'altro
	double a = (ising.Get_temp_high()-ising.Get_temp_low())/(ising.Get_camp());
	for(double temp=ising.Get_temp_high(); temp>=ising.Get_temp_low(); temp-=a){
		double beta = 1.0/temp;
		// ciclo di equilibrazione
		double sigma = 0.02;
	 	ising.Auto_eq(rnd,temp,sigma);
		// calcolo i valori dei 20 blocchi
		cout << "------------Values for temperature = " << temp <<" ----------------------------------" << endl;
		cout << "Energy -- Heat_c -- Magn -- Chi" << endl;
		for(int j=0; j<ising.Get_nblock(); j++){
			ising.Set_u(0);
			ising.Set_u2(0);
			ising.Set_magn(0);
			ising.Set_chi(0);
			for(int i=0; i<ising.Get_nstep(); i++){
				ising.Measure();
				ising.Move(beta,rnd);
			}
			ene[j] = ising.Get_u()/(ising.Get_nstep()*ising.Get_nspin());
			ene2[j] = ising.Get_u2()/(ising.Get_nstep()*ising.Get_nspin());
			heat[j] = pow(beta,2)*(ene2[j]-pow(ene[j],2));
			magnet[j] = ising.Get_magn()/(ising.Get_nstep()*ising.Get_nspin());
			susc[j] = beta*ising.Get_chi()/(ising.Get_nstep()*ising.Get_nspin());
			cout << ene[j] << " " << heat[j] << " " << magnet[j] << " " << susc[j] << endl;
		}
		// calcolo ora errori
		Write << temp << " ";
		calculation(ising,Write,ene,heat,magnet,susc);
	}
	Write.close();
	return 0;
}
