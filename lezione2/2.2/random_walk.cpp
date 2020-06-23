#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include "random.h"
#include "random_walk.h"

using namespace std;

RandomWalk :: RandomWalk(){}

RandomWalk :: ~RandomWalk(){}

void RandomWalk :: SetOrigin(double x, double y, double z){
	X0=x;
	Y0=y;
	Z0=z;
	return;
}

void RandomWalk :: SetLattice(double a){
	A=a;
	return;
}

void RandomWalk :: Walk_cube(Random& rnd, int N, vector<double>& x, vector<double>& y, vector<double>&z){
	int dir; // serve per decidere le il movimento sarà lungo x,y oppure z
	x.push_back(X0);
	y.push_back(Y0);
	z.push_back(Z0); // ho settato l'origine il i v[0] dei vari vettori
	for(int i=1; i<N; i++){
		dir = rnd.Rannyu(1,7); // se dir=1 --> x; dir=2 --> -x ...
		if(dir==1){
			x.push_back(x[i-1]+A);
			y.push_back(y[i-1]);
			z.push_back(z[i-1]);
		}
		if(dir==2){
			x.push_back(x[i-1]-A);
			y.push_back(y[i-1]);
			z.push_back(z[i-1]);
		}
		if(dir==3){
			x.push_back(x[i-1]);
			y.push_back(y[i-1]+A);
			z.push_back(z[i-1]);
		}
		if(dir==4){
			x.push_back(x[i-1]);
			y.push_back(y[i-1]-A);
			z.push_back(z[i-1]);
		}
		if(dir==5){
			x.push_back(x[i-1]);
			y.push_back(y[i-1]);
			z.push_back(z[i-1]+A);
		}
		if(dir==6){
			x.push_back(x[i-1]);
			y.push_back(y[i-1]);
			z.push_back(z[i-1]-A);
		}
	}
	return;
}

void RandomWalk :: Walk_continuum(Random& rnd, int N, vector<double>& x, vector<double>& y, vector<double>& z){
	x.push_back(X0);
	y.push_back(Y0);
	z.push_back(Z0); // ho settato l'origine il i v[0] dei vari vettori
	double theta, phi;
	for(int i=1; i<N; i++){
		theta = rnd.Rannyu(0,M_PI);
		phi = rnd.Rannyu(0,2*M_PI);
		x.push_back(x[i-1]+A*sin(theta)*cos(phi));
		y.push_back(y[i-1]+A*sin(theta)*sin(phi));
		z.push_back(z[i-1]+A*cos(theta));
	}
	return;		
}







