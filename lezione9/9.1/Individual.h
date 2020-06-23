#ifndef __Individual__
#define __Individual__

#include <vector>
#include "random.h"

using namespace std;

class Individual {

private:
  int Number_gene = 32;
  vector<int> Genis;
  vector<double> x, y;
public:
  // constructor
  Individual(Random&);
  // destructor
  ~Individual();
  // setting
  void Set_Number_gene(int a){Number_gene=a;};
  void Set_Gene(int index, int a){Genis[index]=a;};
  // Get_Number_gene
  int Get_Number_gene(){return Number_gene;};
  double Get_gene(int i){return Genis[i];};
  // methods
  double Calc_Fitness_oncirc();
  double Calc_Fitness2_oncirc(vector<double>&, vector<double>&);
};

#endif // __Individual__
