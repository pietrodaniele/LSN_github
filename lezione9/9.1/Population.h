#ifndef __Population__
#define __Population__

#include <vector>
#include "random.h"
#include "Individual.h"
//
using namespace std;

class Population {

private:
  int Number_ind, point_co;
  vector<Individual> individuals;
  double best_fitness;
public:
  // constructor
  Population(Random&, int, vector<double>&, vector<double>&);
  // destructor
  ~Population();
  // settings
  void Set_Number_ind(int a){Number_ind=a;};
  void Set_an_Individual(int a, Individual& Ind);
  void Clean_population(){individuals.clear();};
  void Add_Individual(Individual Ind){individuals.push_back(Ind);};
  // get
  int Get_Number_ind(){return Number_ind;};
  double Get_BestFitnees(vector<double>& x, vector<double>& y){return individuals[Number_ind-1].Calc_Fitness2_oncirc(x,y);};
  Individual Get_an_Individual(int a){return individuals[a];};
  double Average_BHalfpop(vector<double>&, vector<double>&);
  // methods
  int Selection(Random&);
  void New_Generation(vector<Individual>&);
  void Crossover(Random&, Individual&, Individual&);
  void Mutation1(Random&, Individual&); // single inversion
  void Mutation2(Random&, Individual&); // shit i-->i+1
  void Mutation3(Random&, Individual&); // permutation
  void Mutation4(Random&, Individual&); // order inversion
  void In_Fitness_order(vector<double>&, vector<double>&);
  //Population New_Generation(Population*, Population*);
};

#endif // __Population__
