#include "FunzBase.h"
#include <cmath>

class Integranda : public FunzBase {
private:

public:
 	// costruttore e distruttore
	Integranda();
	~Integranda();
 	// operatori
	virtual double Eval(double x) const override {
		return (M_PI/2)*(cos(M_PI*x/2));
	}
	virtual double Eval_sampling(double x) const override {
			return (M_PI/2)*(cos(M_PI*x/2))/(2*(1-x));
	}
};
