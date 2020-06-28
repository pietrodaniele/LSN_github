#include "finance.h"
#include <cmath>

Finance :: Finance(){}

Finance :: ~Finance(){}

double Finance :: europeanCall(double S, double t){
	return S*N(d(S, t, r, sigma, K, T)) - K*exp(-r*(T-t))*N(d(S, t, r, sigma, K, T) -sigma*sqrt(T-t));
};

double Finance :: europeanPut(double S, double t) {
	return S*(N(d(S, t, r, sigma, K, T))-1) - K*exp(-r*(T-t)) * (N(d(S, t, r, sigma, K, T) - sigma*sqrt(T-t))-1);
};

double Finance :: St_direct(double dt, Random& rnd){
	return S0*exp((r-sigma*sigma/2.)*dt + sigma*sqrt(dt)*rnd.Gauss(0,dt));
};

double Finance :: St_discret(double dt, Random& rnd){
	double S=S0, delta_t=double(dt/t_steps), t=0;
	for (int n = 0; n < t_steps; n++) {
		t += delta_t;
		S *= exp((r-sigma*sigma/2.)*delta_t + sigma*sqrt(delta_t)*rnd.Gauss(0,1));
	}
	return S;
};
