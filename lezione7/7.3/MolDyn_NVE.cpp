/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <string>
#include <cstdlib>
#include "MolDyn_NVE.h"
#include "experiment.h"

using namespace std;


int main(int argc,char*argv[]){
  controll(argc);
  if(atoi(argv[1])==0){
  	Input();
  	cout << "First sim" << endl;
  }else{
  	Restart();
  	cout << "Restart" << endl;
  }             //Inizialization
  int nconf = 1;
  int Num_mis =1;
  for(int istep=1; istep <= nstep; ++istep){
     Move();        //Move particles with Verlet algorithm
     if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
     if(istep%10 == 0){
        Measure();     //Properties measurement
        //ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"!
        nconf += 1;
        if(Num_mis==100){
          Print_hist_block(100);
          Clean_hist_block();
          Num_mis=1;
        }else{
          Num_mis++;
        }
     }
  }
  int nmis=nstep/10, L=100, Num=nmis/L; // number of ...
  ConfFinal();         //Write final configuration to restart
  ConfOld();
  Final_g_err(100);
  // this is the part of exercise 4.2+4.3
  string type[4] = {"epot","ekin","temp","etot"};
  Experiment Md;
  double av[Num], av2[Num], sum_prog[Num], sum2_prog[Num], error_prog[Num];
  for(int k=0; k<4; k++){
  	Md.cicleblock(av,av2,type,nmis,L,k);
  	Md.accumulation(sum_prog,av,Num);
  	Md.accumulation(sum2_prog,av2,Num);
  	Md.errorprog(sum_prog,sum2_prog,error_prog,Num);
  	Write(sum_prog,error_prog,type,Num,k,argv[2]);
  }
  return 0;
}

void Write(double* v1, double *v2, string *type, int nmis, int k, char *state){
	ofstream WriteData;
	WriteData.open("ave_"+type[k]+"_"+state+".dat");
	for(int i=0; i<nmis; i++){
		WriteData<<v1[i]<<" "<<v2[i]<<endl;
	}
	WriteData.close();
	return;
};

void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator

  ReadInput.open("input.dat"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables

//Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

//Prepare initial velocities
   cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
   double sumv[3] = {0.0, 0.0, 0.0};
   for (int i=0; i<npart; ++i){
     vx[i] = rand()/double(RAND_MAX) - 0.5;
     vy[i] = rand()/double(RAND_MAX) - 0.5;
     vz[i] = rand()/double(RAND_MAX) - 0.5;

     sumv[0] += vx[i];
     sumv[1] += vy[i];
     sumv[2] += vz[i];
   }
   for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
   double sumv2 = 0.0, fs;
   for (int i=0; i<npart; ++i){
     vx[i] = vx[i] - sumv[0];
     vy[i] = vy[i] - sumv[1];
     vz[i] = vz[i] - sumv[2];

     sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   }
   sumv2 /= (double)npart;

   fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor
   for (int i=0; i<npart; ++i){
     vx[i] *= fs;
     vy[i] *= fs;
     vz[i] *= fs;

     xold[i] = Pbc(x[i] - vx[i] * delta);
     yold[i] = Pbc(y[i] - vy[i] * delta);
     zold[i] = Pbc(z[i] - vz[i] * delta);
   }
   // cleaning hist
   for(int i=0; i<npart; i++){
     hist[i]=0;
   }
   bin_size = (box/2.0)/(double)bins;
   return;
}

void Restart(void){
  cout << "Make sure you've copied first old.final on old.0 and the config.final on config.0" << endl;
  ifstream ReadInput,ReadConf,ReadOld;
  ofstream WritePos;
  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator

  ReadInput.open("input.dat"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables

//Read initial configuration, I need to copy config.final on config.0, but after copying config.0 on old.0
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  ReadOld.open("old.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    ReadOld >> xold[i] >> yold[i] >> zold[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
    xold[i] = xold[i] * box;
    yold[i] = yold[i] * box;
    zold[i] = zold[i] * box;
  }

  ReadConf.close();
  ReadOld.close();
// Now I need to computate r(t+dt) with Verlet integration scheme
  double x_dt[m_part], y_dt[m_part], z_dt[m_part], t=0.0, Temp_dt2, scale_factor, xnew, ynew, znew;
  double fx[m_part], fy[m_part], fz[m_part], vx_dt2[m_part], vy_dt2[m_part], vz_dt2[m_part];
  for(int i=0; i<npart; ++i){ //Force acting on particle i
      fx[i] = Force(i,0);
      fy[i] = Force(i,1);
      fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

      x_dt[i] = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
      y_dt[i] = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
      z_dt[i] = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

      vx_dt2[i] = Pbc(x_dt[i] - x[i])/(delta);
      vy_dt2[i] = Pbc(y_dt[i] - y[i])/(delta);
      vz_dt2[i] = Pbc(z_dt[i] - z[i])/(delta); // v(t+dt/2)

  }
  // Now T(t+dt/2)
  for (int i=0; i<npart; ++i) t += 0.5 * (vx_dt2[i]*vx_dt2[i] + vy_dt2[i]*vy_dt2[i] + vz_dt2[i]*vz_dt2[i]);
  Temp_dt2 = (2.0 / 3.0) * t/(double)npart;
  // By comparing 𝑇(𝑡+𝑑𝑡/2) with the desired/target temperature 𝑇⋆ extract a scaling factor for the velocities and rescale them
  scale_factor =  temp/Temp_dt2;
  cout << endl;
  cout <<"The temperature's scaling factor is = " << scale_factor << endl << endl;
  for(int i=0; i<npart; ++i){
    vx[i] = Pbc(x_dt[i] - xold[i])/(2.0 * delta) * scale_factor;
    vy[i] = Pbc(y_dt[i] - yold[i])/(2.0 * delta) * scale_factor;
    vz[i] = Pbc(z_dt[i] - zold[i])/(2.0 * delta) * scale_factor; // v(t)*scale_factor
    // Ora calcolo con r(t+dt) e v_s(t) r_new(t)
 	xnew = Pbc(x_dt[i]-delta*vx[i]);
    ynew = Pbc(y_dt[i]-delta*vy[i]);
    znew = Pbc(z_dt[i]-delta*vz[i]);
    // E r_new(t)-->r_old; r_dt-->r
    xold[i] = xnew;
    yold[i] = ynew;
    zold[i] = znew;

    x[i] = x_dt[i];
    y[i] = y_dt[i];
    z[i] = z_dt[i];
  }
  // cleaning hist
  for(int i=0; i<npart; i++){
    hist[i]=0;
  }
  bin_size = (box/2.0)/(double)bins;
  return;
}

void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }

  return f;
}

void Measure(){ //Properties measurement
  int bin;
  double v, t, vij, wij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp, Pres;
  double w = 0.0;

  Epot.open("data/output_epot.dat",ios::app);
  Ekin.open("data/output_ekin.dat",ios::app);
  Temp.open("data/output_temp.dat",ios::app);
  Etot.open("data/output_etot.dat",ios::app);
  Pres.open("data/output_pres.dat",ios::app);

  v = 0.0; //reset observables
  t = 0.0;
  // cleaning histogram
  for(int i=0; i<bins; i++) hist[i]=0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
     dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
     dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     if(dr < box/2){
        int m = int(dr/bin_size);
        //cout << m << endl;
        hist[m] += 2;
      }

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
       wij = 1.0/pow(dr,12) - 0.5/pow(dr,6);
//Potential energy
       v += vij;
       w += wij;
     }
    }
    //cout<<1<<endl;
    for(int i=0; i<bins; i++){
      double deltaV =  (4./3.)*M_PI*(pow((i+1)*bin_size,3)-pow((i)*bin_size,3));
      // cout << deltaV << endl;
      hist_block[i] += double(hist[i]/(rho*npart*deltaV));
    }
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);

    stima_pot = v/(double)npart; //Potential energy per particle
    stima_kin = t/(double)npart; //Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total energy per particle
    stima_pres = (48.0 * w / 3.0)/(double)npart; //Pressure per particle

    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;
    Pres << stima_pres << endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();
    Pres.close();

    return;
}

void Print_hist_block(int m){

  ofstream Write;
  Write.open("data/histogram_blocks.0",ios::app);
  for(int f=0; f<bins; f++){
    double r = (f+1)*bin_size;
    //double deltaV =  (4./3.)*pi*(pow(box/2.+r,3)-(pow((box/2.),3)));
    //Write << r << " " << hist_block[f]/((rho*npart)*(deltaV)) << endl;
    Write << r << " " << hist_block[f]/m << endl;
  }
  Write.close();
  return;
};

void Clean_hist_block(){
// cleaning histogram
  for(int i=0; i<bins; i++) hist_block[i]=0;
  return;
};

void Final_g_err(int nblk){

  ifstream Read;
  double sum=0, sum2=0, err=0;
  Read.open("data/histogram_blocks.0");
  double* r = new double[nblk*bins];
  double* H = new double[nblk*bins];
  int m=0;
  while(!Read.eof()){
    Read >> r[m] >> H[m];
    m++;
  }
  cout << 0 << endl;
  ofstream Write("data/final_hist.0");
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
void ConfOld(void){ //Write final configuration
  ofstream WriteOld;

  cout << "Print final configuration to file old.final " << endl << endl;
  WriteOld.open("old.final");

  for (int i=0; i<npart; ++i){
    WriteOld << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  }
  WriteOld.close();
  return;
}

void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

void controll(int a){
	if(a != 3){
		cerr << "./MolDyn_NVE.exe A type_of_state; where A=0/1" << endl;
		exit(3);
	}else{}
	return;
}
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
