#ifndef __Metropolis__
#define __Metropolis__

#include "FunzBase3D.h"
#include "random.h"
#include <iostream>
#include <vector>

using namespace std;

class Metropolis {

private:
	double delta, x_0, y_0, z_0;
protected:

public:
	// constructors
	Metropolis();
  	// destructor
  	~Metropolis();
  	// setting
  	void SetDelta(double);
  	void SetStart(double, double, double);
  	// methods
  	void Algoritm_unif(int n, FunzBase3D&, Random&, vector<double>&, vector<double>&, vector<double>&); // passo il numero di passi da fare, la distribuzione di prob, e la uniform transition probability
									// assumiamo che sia simmetrica
	void Algoritm_gauss(int n, FunzBase3D&, Random&, vector<double>&, vector<double>&, vector<double>&);
	void coord_r(vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&, int);
  	void accumulation(vector<double>&, vector<double>&); // accumulo
  	void errorprog(vector<double>&, vector<double>&, vector<double>&); // calcolo l'errore
};

#endif // __Random__
