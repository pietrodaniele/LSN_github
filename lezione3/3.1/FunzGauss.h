#include "FunzBase.h"
#include <cmath>

class FunzExp : public FunzBase {
private:
	double sigma, mu;
public:
 	// costruttore e distruttore
	FunzGauss();
	~FunzGauss();
	// setto i valori della gaussina
	void Set_mu(double a){mu=a};
	void Set_sigma(double b){sigma=b};
 	// operatori
	virtual double Eval(double x) const override {
		return e^(-(x - mu)^2/(2*(sigma)^2))/(sqrt(2*M_PI)*sigma);
	}
};
