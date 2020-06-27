#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include "random.h"
#include "finance.h"
#include "funzioni.h"
#include "experiment.h"


using namespace std;

// objects
Experiment lab;
Finance market;
Random rnd;
// parameters
double S0=100, T=1, K=100, r=0.1, sigma=0.25;
unsigned int t_steps=100, M=100000, n=1000, l=M/n;
vector<double> EC_direct, EC_discret, EP_direct, EP_discret;

int main (){
	cout << "MC steps M = "<< M << endl;
	cout << "# of blocks N = " << n << endl;
	cout << "# of MC steps in each block = " << l << endl;
	cout << endl;
	cout << "---------------------------------------" << endl;
	cout << endl;
	cout << "Market's parameters:" << endl;
	cout << "		- asset price at t=0: S(0)=" << S0 << endl;
	cout << "		- delivery time: T=" << T << endl;
	cout << "		- strike price: K="<< K << endl;
	cout << "		- risk-free interest rate: r=" << r << endl;
	cout << "		- volatility: sigma=" << sigma << endl;
	cout << endl;
	// settings
	market.SetS0(S0);
  market.SetT(T);
  market.SetK(K);
	market.Setr(r);
	market.SetSigma(sigma);
	market.Set_t_steps(t_steps);
	// creating data
	for(unsigned int d=0; d<M; d++){
		EC_direct.push_back(market.europeanCall(market.St_direct(),T));
		EC_discret.push_back(market.europeanCall(market.St_discret(),T));
		EP_direct.push_back(market.europeanPut(market.St_direct(),T));
		EP_discret.push_back(market.europeanPut(market.St_discret(),T));
	}
	lab.Block_prog_ave_print(EC_direct,"EC_direct",n,l);
	cout << "European Call with S(t) direct: done" << endl;
	lab.Block_prog_ave_print(EC_discret,"EC_discret",n,l);
	cout << "European Call with S(t) discret: done" << endl;
	lab.Block_prog_ave_print(EP_direct,"EP_direct",n,l);
	cout << "European Put with S(t) direct: done" << endl;
	lab.Block_prog_ave_print(EP_discret,"EP_discret",n,l);
	cout << "European Put with S(t) discret: done" << endl;
	// clean
	EC_direct.clear();
	EC_discret.clear();
	EP_direct.clear();
	EP_discret.clear();
	return 0;
}
