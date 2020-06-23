#include "FunzBase.h"
#include <cmath>

class FunzLnz : public FunzBase {
private:
	double T=1;
public:
 	// costruttore e distruttore
	FunzLnz();
	~FunzLnz();
 	// operatori
	virtual double Eval(double x) const override {
		return T*tan(M_PI*(x-1/2));
	}
};
