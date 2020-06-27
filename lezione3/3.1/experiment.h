#ifndef __Experiment__
#define __Experiment__

#include <iostream>
using namespace std;

class Experiment{

private:
  vector<double> av, av2, sum_prog, sum2_prog, error_prog;
protected:

public:
  // constructors
  Experiment();
  // destructor
  ~Experiment();
  // methods: costruisco i vari metodi. Uno per fare il vettore sommato. Uno per gli accumuli ecc...
  void cicleblock(vector<double>&, int, int);
  void accumulation(vector<double>&, vector<double>&);
  void errorprog();
  void Block_prog_ave_print(vector<double>&, string, int, int);
};

#endif // __Experiment__
