#ifndef __Metropolis__
#define __Metropolis__

#include <vector>
#include "FunzBase.h"
#include "random.h"

using namespace std;

class Metropolis {

private:
	double delta, x_0;
protected:

public:
	// constructors
	Metropolis();
  // destructor
  ~Metropolis();
  // setting
  void SetDelta(double);
  void SetStart(double);
  // methods
  void Algoritm_unif(int, FunzBase&, Random&, vector<double>&);// passo il numero di passi da fare, la distribuzione di prob, e la uniform transition probability
	void Algoritm_unif_Try(int, FunzBase&, Random&, vector<double>&);
	// integral with data blocking
	void cicleblock_integral(vector<double>&, vector<double>&, int, int, FunzBase&);
};

#endif // __Random__


// includo il random walk per il T(x'|x_n); devo però implementare il fatto che sia continuo qui comprende tutti i punti all'interno
// del delta. Al di fuori di delta sarà nulla
