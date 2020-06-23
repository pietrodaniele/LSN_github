#include "FunzBase.h"
#include <cmath>

class FunzExp : public FunzBase {
private:
	double lambda=1;
public:
 	// costruttore e distruttore
	FunzExp();
	~FunzExp();
 	// operatori
	virtual double Eval(double x) const override {
		return (-1/lambda)*log(1-x);
	}
};
