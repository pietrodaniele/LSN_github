#ifndef __Experiment__
#define __Experiment__

#include <iostream>
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
  void cicleblock(vector<double>&, vector<double>&, Random, int, int);
  void accumulation(vector<double>&, vector<double>&, int);
  void errorprog(vector<double>&, vector<double>&, vector<double>&, int);
};

#endif // __Experiment__

