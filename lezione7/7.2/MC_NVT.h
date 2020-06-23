#ifndef __MC_NVT__
#define __MC_NVT__

#include <vector>
#include "random.h"

class Mc_nvt {

private:
  double temp, beta, rho, rcut, delta, bin_size, box, vol, vtail, ptail, Uint, Pres, Virial, dR;
  const double pi=3.1415927;
  int npart, nblk, nstep, bins=100, number_s, eq_gas=1500, eq_liquid=1000, eq_solid=750, accepted, attempted;
  int seed[4];
  double* hist_final = new double[bins];
  double* hist_block = new double[bins];
  double* hist = new double[bins];
  double* x = new double[108];
  double* y = new double[108];
  double* z = new double[108];
protected:

public:
  // constructors
	Mc_nvt();
  // destructor
  ~Mc_nvt();
  // setting

  // get
  int Get_nblk(){return nblk;};
  int Get_nstep(){return nstep;};
  int Get_acc(){return accepted;};
  int Get_att(){return attempted;};
  // methods
  void Input(Random&, int);
  void Move(Random&);
  void Measure();
  void Equilibration(Random&, int);
  void Print_hist_block();
  void Clean_hist_block();
  void Final_g_err();
  double Pbc(double);
  double Boltzmann(double, double, double, int);
};

#endif // __MC_NVT__
