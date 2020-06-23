#include "FunzBase.h"

class FunzSigma : public FunzBase {
private:

public:
 	// costruttore e distruttore
	FunzSigma();
	~FunzSigma();
 	// operatori
	virtual double Eval(double x) const override {
		return (x-0.5)*(x-0.5);
	}
};
