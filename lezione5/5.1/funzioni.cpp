#include "funzioni.h"
#include "random.h"
#include <string>

void createrandom(Random& rnd){
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;
   rnd.SaveSeed();
};

void accettanza_uniform(double alfa, vector<double>& x, vector<double>& y, vector<double>& z, Random& rnd, double& acc,\
						 int i, double x_try, double y_try, double z_try){
	double r=rnd.Rannyu();
	if(r<=alfa){
		x.push_back(x_try);
		y.push_back(y_try);
		z.push_back(z_try);
		acc++;
	}else{
		x.push_back(x[i-1]);
		y.push_back(y[i-1]);
		z.push_back(z[i-1]);
	}
	return;
};

void write_100(vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& sum, vector<double>& error){
	ofstream write;
	write.open("data/hydro_100.dat");
	for(unsigned int h=0; h<z.size(); h++){
		write << x[h] << " "<< y[h] << " "<< z[h]  << endl;
	}
	write.close();
	write.open("data/r_100.dat");
	for(unsigned int h=0; h<sum.size(); h++){
		write << h << " " << sum[h] << " "<< error[h]  << endl;
	}
	write.close();
	return;
}

void write_100gauss(vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& sum, vector<double>& error){
	ofstream write;
	write.open("data/hydro_100gauss.dat");
	for(unsigned int h=0; h<z.size(); h++){
		write << x[h] << " "<< y[h] << " "<< z[h]  << endl;
	}
	write.close();
	write.open("data/r_100gauss.dat");
	for(unsigned int h=0; h<sum.size(); h++){
		write << h << " " << sum[h] << " "<< error[h]  << endl;
	}
	write.close();
	return;
}

void write_210(vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& sum, vector<double>& error){
	ofstream write;
	write.open("data/hydro_210.dat");
	for(unsigned int h=0; h<z.size(); h++){
		write << x[h] << " "<< y[h] << " "<< z[h]  << endl;
	}
	write.close();
	write.open("data/r_210.dat");
	for(unsigned int h=0; h<sum.size(); h++){
		write << h << " " << sum[h] << " "<< error[h]  << endl;
	}
	write.close();
}

void write_210gauss(vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& sum, vector<double>& error){
	ofstream write;
	write.open("data/hydro_210gauss.dat");
	for(unsigned int h=0; h<z.size(); h++){
		write << x[h] << " "<< y[h] << " "<< z[h]  << endl;
	}
	write.close();
	write.open("data/r_210gauss.dat");
	for(unsigned int h=0; h<sum.size(); h++){
		write << h << " " << sum[h] << " "<< error[h]  << endl;
	}
	write.close();
	return;
}

void clean(vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& v4, vector<double>& v5, vector<double>& v6, vector<double>& v7, vector<double>& v8){
	x.clear();
	y.clear();
	z.clear();
	v4.clear();
	v5.clear();
	v6.clear();
	v7.clear();
	v8.clear();
	return;
};

double minimum(double a, double b){
	if(a<=b){
		return a;
	}else{
		return b;
	}
};
