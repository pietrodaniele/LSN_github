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
  void cicleblock(vector<double>&, Random&, int, int); // ho tolto il av2 !!!!!!!!!!!!!!
  void cicleblock_funz(vector<double>&, Random&, int, int, FunzBase& );// pure qui!!!!!!!!!!!!
};

#endif // __Experiment__

