#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include "random.h"
#include "experiment.h"

using namespace std;

Experiment :: Experiment(){}

Experiment :: ~Experiment(){}

void Experiment :: cicleblock(vector<double>& av, vector<double>& av2, Random& rnd, int N, int L){
	vector<double> r;
	double sum;
	int k=0;
	for(int i=0; i<N*L; i++){
		r.push_back(rnd.Rannyu());
	}
	for(int i=0; i<N; i++){
		sum=0;
		for(int j=0; j<L; j++){
			k=j+i*L;
			sum += r[k];
		}
		av.push_back(sum/L);
		av2.push_back(av[i]*av[i]);
	}
	return;
}

void Experiment :: accumulation(vector<double>& sum_prog, vector<double>& av, int N){
	for(int i=0; i<N; i++){
		sum_prog.push_back(0);
	}
	for(int i=0; i<N; i++){
		for(int j=0; j<i+1; j++){
			sum_prog[i]+=av[j];
		}
		sum_prog[i]=sum_prog[i]/(i+1);
	}
	return;	
}

void Experiment :: errorprog(vector<double>& sum_prog, vector<double>& sum2_prog, vector<double>& error_prog, int N){
	for(int i=0; i<N; i++){
		if(i==0){
			error_prog.push_back(0);
		 }else{
			error_prog.push_back(sqrt((sum2_prog[i]-(sum_prog[i]*sum_prog[i]))/i));
		}
	}
	return;
}
	
void Experiment :: cicleblock_funz(vector<double>& av, vector<double>& av2, Random& rnd, int N, int L, FunzBase& f){
		vector<double> r;
	double sum;
	int k=0;
	for(int i=0; i<N*L; i++){
		r.push_back(rnd.Rannyu());
	}
	for(int i=0; i<N; i++){
		sum=0;
		for(int j=0; j<L; j++){
			k=j+i*L;
			sum += f.Eval(r[k]);
		}
		av.push_back(sum/L);
		av2.push_back(av[i]*av[i]);
	}
	return;
}

void Experiment :: chiquadro(Random& rnd, int M, int n, vector<double>& chiquadro){
	double sum, c, bin_sx, bin_dx;
	int b;
	int* n_i = new int[M]; // vector of the frequency
	for(int i=0; i<M; i++){
		n_i[i]=0;
	}
	// the following "for" calculate M chiquadro's values
	for(int k=0; k<M; k++){
		sum=0;
		// divide [0,1] in M bins. Now I need to calculate the frequency n_i
		for(int j=0; j<n; j++){
			c=rnd.Rannyu();
			for(double i=0; i<M; i++){
				bin_sx=(i/M);
				bin_dx=((i+1)/M);
				//cout<<bin_sx<<" < "<<bin_dx<<endl;
				if(c>=bin_sx && c<bin_dx){
					b=i;
					n_i[b]+=1;
				}else{}
			}
		}
		// now I need to calculate the chiquadro_k, so I write a "for" that make the sum
		for(int h=0; h<M; h++){
			sum+=(pow(n_i[h]-(n/M),2)/(n/M));
		}
		chiquadro.push_back(sum);
		for(int i=0; i<M; i++){
				n_i[i]=0;
		}
	}

	return;
}





