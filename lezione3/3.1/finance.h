#ifndef __Finance__
#define __Finance__

class Finance {

private:
	double S0=100, T=1, K=100, r=0.1, sigma=0.25;

protected:

public:
	// constructors
	Finance();
    // destructor
    ~Finance();
    /*// setto i parametri
    void SetS0(double A){S0=A;};
    void SetT(double A){T=A;};
    void SetK(double K){S0=S};*/// aggiungere i get per i private in modo da stampare i parametri
    // methods
	double europeanCall();
	double europeanPut();
	double St_unif(double , double , double );
	double St_disc();
	
};

#endif // __Finance__


// primo modo simuliamo il prezzo finale direttamente
// secondo modo dividiamo 
