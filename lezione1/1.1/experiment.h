#ifndef __Experiment__
#define __Experiment__

#include <iostream>
#include "FunzBase.h"

using namespace std;

class Experiment{

private:

protected:

public:
  // constructors
  Experiment();
  // destructor
  ~Experiment();
  // methods: costruisco i vari metodi. Uno per fare il vettore sommato. Uno per gli accumuli ecc...
  void cicleblock(vector<double>&, vector<double>&, Random&, int, int);
  void cicleblock_funz(vector<double>&, vector<double>&, Random&, int, int, FunzBase& );
  void accumulation(vector<double>&, vector<double>&, int);
  void errorprog(vector<double>&, vector<double>&, vector<double>&, int);
  void chiquadro(Random&, int, int , vector<double>&);
};

#endif // __Experiment__

