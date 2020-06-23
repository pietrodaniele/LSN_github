#ifndef __RandomWalk__
#define __RandomWalk__

#include <iostream>
#include <vector>

using namespace std;

class RandomWalk{ // è un random walk 3D

private:
	double X0, Y0, Z0, A; // metto in private le coordinate dell'origine e poi ne metodi aggiungo il settaggio
protected:

public:
  	// constructors
	RandomWalk();
  	// destructor
  	~RandomWalk();
  	// methods
  	void SetOrigin(double,double,double);
  	void SetLattice(double);
  	void Walk_cube(Random&, int, vector<double>&, vector<double>&, vector<double>&); // random per dare la direzione, int per il numero di passi
  	void Walk_continuum(Random&, int, vector<double>&, vector<double>&, vector<double>&); 
};

#endif // __RandomWalk__


