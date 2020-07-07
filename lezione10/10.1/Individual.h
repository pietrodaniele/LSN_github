#ifndef __Individual__
#define __Individual__

#include <vector>
#include "random.h"

using namespace std;

class Individual {

private:
  int Number_gene;
  vector<int> Genis;
public:
  // constructor
  Individual(Random&, int);
  // destructor
  ~Individual();
  // setting
  void Set_Number_gene(int a){Number_gene=a;};
  void Set_Gene(int index, int a){Genis[index]=a;};
  void Set_Genis(vector<int>);
  // Get_Number_gene
  int Get_Number_gene(){return Number_gene;};
  double Get_gene(int i){return Genis[i];};
  vector<int> Get_Genis(){return Genis;};
  // methods
  double Calc_Fitness2(vector<double>&, vector<double>&);
  void Evolve(Random&);
  void Write_path(int, vector<double>&, vector<double>&);
  // mutation
  void Mutation1(int, int);
  void Mutation2(int, int);
  void Mutation3(int, int);
  void Mutation4(int, int);
  void Print();
};

#endif // __Individual__
