#include "FunzBase.h"
#include <cmath>

class Funz_Try : public FunzBase {
private:
	double mu, sigma;
public:
	// settings
	void Set_Mu(double x){mu=x;};
	void Set_Sigma(double x){sigma=x;};
	// get
	double Get_Mu(){return mu;};
	double Get_Sigma(){return sigma;};
 	// costruttore e distruttore
	Funz_Try();
	~Funz_Try();
 	// operatori
	virtual double Eval(double x) const override {
		//double norm = //2*pow(2*M_PI,0.5)/pow(1/(sigma*sigma),0.5);
		//double norm = 0.5*sqrt(M_PI)*sigma*((2*exp(-(mu*mu)/(sigma*sigma))*erf(x/sigma) + erf((x - mu)/sigma) + erf((x + mu)/sigma)));
		return (exp(-((x-mu)*(x-mu))/(2*sigma*sigma))+exp(-((x+mu)*(x+mu))/(2*sigma*sigma)));
	};
	virtual double Eval2(double x) const override {
		//double norm = //2*pow(2*M_PI,0.5)/pow(1/(sigma*sigma),0.5);
		//double norm = 0.5*sqrt(M_PI)*sigma*((2*exp(-(mu*mu)/(sigma*sigma))*erf(x/sigma) + erf((x - mu)/sigma) + erf((x + mu)/sigma)));
		return pow((exp(-((x-mu)*(x-mu))/(2*sigma*sigma))+exp(-((x+mu)*(x+mu))/(2*sigma*sigma))),2);
	};
	virtual double EvalD2(double x) const override {
		return -(exp(-(x-mu)*(x-mu)/(2*sigma*sigma)))/(sigma*sigma) \
						-(exp(-(x+mu)*(x+mu)/(2*sigma*sigma)))/(sigma*sigma) \
						+(pow(x-mu,2))*(exp(-(x-mu)*(x-mu)/(2*sigma*sigma)))/(pow(sigma,4)) \
						+(pow(x+mu,2))*(exp(-(x+mu)*(x+mu)/(2*sigma*sigma)))/(pow(sigma,4));
	};
};
