#ifndef __SimAnneling__
#define __SimAnneling__

#include <vector>
#include "random.h"
#include "Individual.h"

using namespace std;

class SimAnneling {

private:
  double start_Temp, deltaTemp, absTemp, Temp;
public:
  // constructor
  SimAnneling();
  // destructor
  ~SimAnneling();
  // setting
  void Set_start_Temp(double T){start_Temp=T;};
  void Set_dT(double dt){deltaTemp=dt;};
  void Set_absTemp(double abTEmp){abTEmp=absTemp;};
  // Get_Number_gene
  double Get_Temp(){return start_Temp;};
  double Get_dT(){return deltaTemp;};
  double Get_absTemp(){return absTemp;};
  // methods
  void Anneal(Individual&, Random&, vector<double>&, vector<double>&, int);
};

#endif // __SimAnneling__
