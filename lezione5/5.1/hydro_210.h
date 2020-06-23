#include "FunzBase3D.h"
#include <cmath>

class Hydro_210 : public FunzBase3D {
private:

public:
 	// costruttore e distruttore
	Hydro_210();
	~Hydro_210();
 	// operatori
	virtual double Eval(double x, double y, double z) const override {
		double a_0 = 0.0529; 
		double r = pow(x*x + y*y + z*z, 0.5);
		return pow(((pow(a_0,(-5.0/2.0)))/8)*(pow(2.0/M_PI,0.5))*r*exp(-r/(2*a_0))*z/r,2);
	}
};
