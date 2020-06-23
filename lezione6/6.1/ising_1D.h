#ifndef __Ising_1D__
#define __Ising_1D__

#include <vector>
#include "random.h"

class Ising_1D {

private:
	double temp_low, beta_low, temp_high, beta_high, J, h, metro;
	double u, u2, magn, chi; // valori instantanei
	int nspin, nblock, nstep, camp;
	double* s = new double[50];
protected:

public:
	// constructors
	Ising_1D();
  	// destructor
  	~Ising_1D();
  	// setting
  	void Input(Random&);
  	void Set_u(double e){u=e;};
  	void Set_u2(double e2){u2=e2;};
  	void Set_magn(double m){magn=m;};
  	void Set_chi(double m){chi=m;};
  	// get
  	int Get_nstep(){return nstep;};
  	int Get_nspin(){return nspin;};
  	int Get_nblock(){return nblock;};
  	int Get_camp(){return camp;};
  	double Get_temp_low(){return temp_low;};
  	double Get_temp_high(){return temp_high;};
  	double Get_metro(){return metro;};
  	double Get_J(){return J;};
  	double Get_h(){return h;};
  	double Get_u(){return u;};
  	double Get_u2(){return u2;};
  	double Get_magn(){return magn;};
  	double Get_chi(){return chi;};
  	// methods
  	void Move(double, Random&);
  	void Measure();
  	void Auto_eq(Random&, double, double);
  	void av2(double*, double*);
  	void accumulation(double*, double*,double*, double*, double*); // accumulo e calcolo errore
};

#endif // __Random__
