#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include "metropolis.h"

using namespace std;

int main() {


	Input();
	

	double x = 0.; 

	double integral = 0.;
	int attempted = 0;
	int accepted = 0;


	double xNew = 0.;
	double A = 0.;

	for (int i=0.; i<nstep; i++) {
		xNew = x + T(x);
		attempted ++;

		if(p(x)==0.) A = 1.;
		else A = min(1.,p(xNew)/p(x));

		if (rnd.Rannyu() <= A){
			x = xNew;
			integral += x;
			accepted++;
		}
	
	}
	cout << "integrale di x: " <<  integral/accepted << endl;
	cout << "rate accettazione: " << double(accepted)/double(attempted) << endl;


	rnd.SaveSeed();
	return 0;
}



void Input(void) {
	
	// Read seed for random numbers
	int p1, p2;
	ifstream Primes("Primes");
	Primes >> p1 >> p2;
	Primes.close();

	ifstream input("seed.in");
	input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
	rnd.SetRandom( seed, p1, p2 );
	input.close();
	

}


double p(double x) {

	if(x>0.) return exp(-2*x);
	else return 0.;

}


double T(double x) {
	return rnd.Rannyu() * delta - delta/2.;
}
