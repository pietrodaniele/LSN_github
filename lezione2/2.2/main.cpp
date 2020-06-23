#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "random.h"
#include "random_walk.h"
#include "funzioni.h"

using namespace std;

Random rnd;
RandomWalk drunk;
vector<double> x, y, z, r, error_r; // su r andrò a salvare 
int N=100, M=10000; // numero di random step e numero di random walk
double x_0=0, y_0=0, z_0=0, a=1; //informazioni sull'origine e sul lattice
 
int main (){
	// setto il generatore del random, il mio cube lattice
    createrandom(rnd);
    drunk.SetOrigin(x_0,y_0,z_0);
    drunk.SetLattice(a);
    // calcolo <|r^2|> e le sue incertezze per il cammino sul reticolo 
    r2_mean_lattice(N,M,rnd,drunk,r,x,y,z);
    error(r,error_r);
    // scrivo il file e poi pulisco i miei vettori per riusarli
    write_lattice(r,error_r);
    clean3(x,y,z);
    clean2(r,error_r);
    // calcolo <|r^2|> e le sue incertezze per il cammino nel continuo
    r2_mean_continuum(N,M,rnd,drunk,r,x,y,z);
    error(r,error_r);
    // scrivo il file e poi pulisco i miei vettori
    write_continuum(r,error_r);
    clean3(x,y,z);
return 0;
}
