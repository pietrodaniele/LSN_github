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
  void cicleblock(double *, double *, string *, int , int, int); // legge il file e calcola av e av2
  void accumulation(double *, double *, int); // accumulo qui per non appesantire troppo cicleblock
  void errorprog(double*, double *, double *, int); // calcolo l'errore
};

#endif // __Experiment__


