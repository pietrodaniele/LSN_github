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

void Experiment :: cicleblock(vector<double>& r, int n, int l){

	double sum;
	int k=0;
	for(int i=0; i<n; i++){
		sum=0;
		for(int j=0; j<l; j++){
			k=j+i*l;
			sum += r[k];
			//cout << k << endl;
		}
		av.push_back(sum/l);
		av2.push_back(av[i]*av[i]);
	}
	return;
}

void Experiment :: accumulation(vector<double>& a, vector<double>& b){
	for(unsigned int i=0; i<b.size(); i++){
		a.push_back(0);
	}
	for(unsigned int i=0; i<b.size(); i++){
		for(unsigned int j=0; j<i+1; j++){
			a[i]+=b[j];
		}
		a[i]=a[i]/(i+1);
	}
	return;
}

void Experiment :: errorprog(){
	for(unsigned int i=0; i<sum_prog.size(); i++){
		if(i==0){
			error_prog.push_back(0);
		 }else{
			error_prog.push_back(sqrt((sum2_prog[i]-(sum_prog[i]*sum_prog[i]))/i));
		}
	}
	return;
}

void Experiment :: Block_prog_ave_print(vector<double>& r, string filename, int n, int l){
	ofstream Write;
	Write.open("data/"+filename+".dat");
	cicleblock(r,n,l);
	accumulation(sum_prog,av);
	accumulation(sum2_prog,av2);
	errorprog();
	for(unsigned int i=0; i<error_prog.size(); i++){
		Write<<i<<" "<<sum_prog[i]<<" "<<error_prog[i]<<endl;
	}
	Write.close();
	av.clear();
	av2.clear();
	sum_prog.clear();
	sum2_prog.clear();
	error_prog.clear();
}
