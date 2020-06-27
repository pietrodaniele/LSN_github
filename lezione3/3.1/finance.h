#ifndef __Finance__
#define __Finance__

#include "funzioni.h"

class Finance {

private:
	/*double S0=100, T=1, K=100, r=0.1, sigma=0.25;
	int t_steps=100;*/
	double S0, T, K, r, sigma;
	int t_steps;
protected:

public:
	// constructors
	Finance();
  // destructor
  ~Finance();
  // setto i parametri
  void SetS0(double A){S0=A;};
  void SetT(double A){T=A;};
  void SetK(double A){K=A;};
	void Setr(double A){r=0.1;};
	void SetSigma(double A){sigma=0.25;};
	void Set_t_steps(double A){t_steps=100;};
  // methods finance
	double europeanCall(double, double);
	double europeanPut(double, double);
	double St_direct();
	double St_discret();
};

#endif // __Finance__


// primo modo simuliamo il prezzo finale direttamente
// secondo modo dividiamo
