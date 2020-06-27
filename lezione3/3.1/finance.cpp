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

double Finance :: St_direct(){
	return 1;
};

double Finance :: St_discret(){
	return 1;
};
