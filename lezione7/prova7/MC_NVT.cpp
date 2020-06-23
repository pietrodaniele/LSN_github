#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "funzioni.h"
#include "random.h"
#include "MC_NVT.h"

using namespace std;

Mc_nvt :: Mc_nvt(){}

Mc_nvt :: ~Mc_nvt(){}

void Mc_nvt :: Input(Random& rnd, int type){
  ifstream ReadInput,ReadConf;
  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T " << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;
  //Read seed for random numbers
  int p1, p2;
  ifstream Primes("Primes");
  Primes >> p1 >> p2 ;
  Primes.close();

  ifstream input("seed.in");
  input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed,p1,p2);
  input.close();
  //Read input informations
  number_s = type;
  if(number_s==1){
    cout << "-------------------GAS-------------------" << endl;
    ReadInput.open("input.gas");
  }if(number_s==2){
    cout << "-------------------LIQUID-------------------" << endl;
    ReadInput.open("input.liquid");
  }if(number_s==3){
    cout << "-------------------SOLID-------------------" << endl;
    ReadInput.open("input.solid");
  }/*else{
      cerr << "Wrong number of the state (1,2,3)" << endl;
      exit(2);
  }*/

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  box = pow(vol,1.0/3.0);
  cout << "Volume of the simulation box = " << vol << endl;
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  cout << "Cutoff of the interatomic potential = " << rcut << endl << endl;

  //Tail corrections for potential energy and pressure
  vtail = (8.0*pi*rho)/(9.0*pow(rcut,9)) - (8.0*pi*rho)/(3.0*pow(rcut,3));
  ptail = (32.0*pi*rho)/(9.0*pow(rcut,9)) - (16.0*pi*rho)/(3.0*pow(rcut,3));
  cout << "Tail correction for the potential energy = " << vtail << endl;
  cout << "Tail correction for the virial           = " << ptail << endl;

  ReadInput >> delta;

  ReadInput >> nblk;

  ReadInput >> nstep;

  cout << "The program perform Metropolis moves with uniform translations" << endl;
  cout << "Moves parameter = " << delta << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();

  //Prepare arrays for measurements
  /*iv = 0; //Potential energy
  iw = 1; //Virial

  n_props = 2; //Number of observables

  //measurement of g(r)
  nbins = 100;*/
  bin_size = (box/2.0)/(double)bins;
  //Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = Pbc( x[i] * box );
    y[i] = Pbc( y[i] * box );
    z[i] = Pbc( z[i] * box );
  }
  ReadConf.close();
  // Equilibration
  Equilibration(rnd,type);
  //Evaluate potential energy and virial of the initial configuration
  Measure();
  // accepted
  accepted = 0;
  attempted = 0;
  //Print initial values for the potential energy and virial
  cout << "Initial potential energy (after eq)(with tail corrections) = " << Uint << endl;
  cout << "Virial                   (after eq)(with tail corrections) = " << Virial << endl;
  cout << "Pressure                 (after eq)(with tail corrections) = " << Pres << endl << endl;
  return;
}

void Mc_nvt :: Move(Random& rnd){
  int o;
  double p, energy_old, energy_new;
  double xold, yold, zold, xnew, ynew, znew;

  for(int i=0; i<npart; ++i){
  //Select randomly a particle (for C++ syntax, 0 <= o <= npart-1)
    o = (int)(rnd.Rannyu()*npart);
  //Old
    xold = x[o];
    yold = y[o];
    zold = z[o];
    energy_old = Boltzmann(xold,yold,zold,o);
  //New
    xnew = Pbc( x[o] + delta*(rnd.Rannyu() - 0.5) );
    ynew = Pbc( y[o] + delta*(rnd.Rannyu() - 0.5) );
    znew = Pbc( z[o] + delta*(rnd.Rannyu() - 0.5) );
    energy_new = Boltzmann(xnew,ynew,znew,o);
  //Metropolis test
    p = exp(beta*(energy_old-energy_new));
    if(p >= rnd.Rannyu()){
    //Update
       x[o] = xnew;
       y[o] = ynew;
       z[o] = znew;

       accepted = accepted + 1.0;
    }
    attempted = attempted + 1.0;
  }
};

void Mc_nvt :: Measure(){
  double v = 0.0, w = 0.0;
  double vij, wij;
  double dx, dy, dz, dr;
// cleaning histogram
// for(int i=0; i<bins; i++) hist[i]=0;
//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){
// distance i-j in pbc
     dx = Pbc(x[i] - x[j]);
     dy = Pbc(y[i] - y[j]);
     dz = Pbc(z[i] - z[j]);

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);
//update of the histogram of g(r)
    if(dr < box/2){
      hist_block[int(dr/bin_size)] += 2;
    }


     if(dr < rcut)
     {
       vij = 1.0/pow(dr,12) - 1.0/pow(dr,6);
       wij = 1.0/pow(dr,12) - 0.5/pow(dr,6);

// contribution to energy and virial
       v += vij;
       w += wij;
     }
    }
  }
  /*ofstream WriteInsta;
  string type[] = {"gas","liquid","solid"};
  WriteInsta.open("insta_values_"+type[number_s-1]+".0",ios::app);
  Uint = 4.0 * v/(double)npart+vtail;
  Virial = ((48.0 * w / 3.0)/(double)npart) + ptail;
  Pres = rho *  temp + ((48.0 * w / 3.0) + (double)npart * ptail) / vol;
  WriteInsta << Uint << " " <<  Pres << endl;
  WriteInsta.close();*/
  ofstream WriteInsta;
  string type[] = {"gas","liquid","solid"};
  WriteInsta.open("data/insta_values_"+type[number_s-1]+".0",ios::app);
  double app_U = 4.0 * v;
  double app_V = 48.0 * w / 3.0;
  Uint = app_U/(double)npart+vtail;
  Virial = app_V/(double)npart + ptail;
  Pres = rho * temp + (app_V + (double)npart * ptail) / vol;
  WriteInsta << Uint << " " << Pres << endl;
  WriteInsta.close();
};

void Mc_nvt :: Equilibration(Random& rnd, int state_of_the_system){
  // reaching the equilibrium
  if(state_of_the_system==1){
    for(int h=0; h<eq_gas; h++){
      Move(rnd);
    }
    cout << "Gas' equilibrium reached up after " << eq_gas << " MC steps" << endl;
  }if(state_of_the_system==2){
    for(int h=0; h<eq_liquid; h++){
      Move(rnd);
    }
    cout << "Liquid's equilibrium reached up after " << eq_liquid << " MC steps" << endl;
  }if(state_of_the_system==3){
    for(int h=0; h<eq_solid; h++){
      Move(rnd);
    }
    cout << "Solid's equilibrium reached up after " << eq_solid << " MC steps" << endl;
  }
  return;
};

void Mc_nvt :: Print_hist_block(){
  accepted = 0;
  attempted = 0;
  ofstream Write;
  string type[] = {"gas","liquid","solid"};
  Write.open("data/histogram_"+type[number_s-1]+".0",ios::app);
  for(double r=bin_size; r<=(box/2.); r+=bin_size){
    double deltaV =  (4./3.)*pi*((pow(r+(box/2.),3))-pow(box/2.,3));
    Write << r << " " << hist_block[int(r/bin_size)]/((rho*npart*nstep)*(deltaV)) << endl;
  }
  return;
};

void Mc_nvt :: Clean_hist_block(){
// cleaning histogram
  for(int i=0; i<bins; i++) hist_block[i]=0;
  return;
};

void Mc_nvt :: Final_g_err(){
  string type[] = {"gas","liquid","solid"};
  ifstream Read;
  double sum, sum2, err;
  Read.open("data/histogram_"+type[number_s-1]+".0");
  double* r = new double[nblk*nstep];
  double* H = new double[nblk*nstep];
  int m=0;
  while(!Read.eof()){
    Read >> r[m] >> H[m];
    m++;
  }
  ofstream Write("data/final_hist_"+type[number_s-1]+".0");
  for(int i=0; i<m/(double)nblk -1; i++){
    for(int j=0; j<nblk; j++){
      sum += H[i+m/nblk];
      sum2 += pow(H[i+m/nblk],2);
    }
    sum /= nblk;
    sum2 /= nblk;
    err=sqrt(abs((sum2 - sum*sum))/(double)nblk);
    Write << r[i] << " " << sum << " " << err << endl;
    sum=0;
    sum2=0;
    err=0;
  }
  Read.close();
  Write.close();
  return;
};



double Mc_nvt :: Pbc(double r){
  return r - box * rint(r/box);
};

double Mc_nvt :: Boltzmann(double xx, double yy, double zz, int ip)
{
  double ene=0.0;
  double dx, dy, dz, dr;

  for (int i=0; i<npart; ++i)
  {
    if(i != ip)
    {
// distance ip-i in pbc
      dx = Pbc(xx - x[i]);
      dy = Pbc(yy - y[i]);
      dz = Pbc(zz - z[i]);

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);

      if(dr < rcut)
      {
        ene += 1.0/pow(dr,12) - 1.0/pow(dr,6);
      }
    }
  }

  return 4.0*ene;
};
