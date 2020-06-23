#ifndef __Experiment__
#define __Experiment__

#include <iostream>
#include <vector>
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
  void cicleblock_integral(vector<double>&, vector<double>&, Random&, int, int, int, int, FunzBase& ); // è stato modificato!!! rispetto al solito
  void valoriquad(vector<double>&, vector<double>&); // valori quad qui per non appesantire troppo cicleblock_funz
  void accumulation(vector<double>&, vector<double>&); // accumulo qui per non appesantire troppo cicleblock_funz
  void errorprog(vector<double>&, vector<double>&, vector<double>&); // calcolo l'errore
};

#endif // __Experiment__


