#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <armadillo>
#include "../random/random.h"
#include <functional> 
#include <iomanip> 

using namespace arma;
using namespace std;

int seed[4];
Random rnd;

// Simulation
int nstep=3000000;
double delta = 1.6;

double p(double,double,double);
double T(double);

int main() {

   int seed=23;
   cout << "seme: "<<endl;
   
   rnd.SetSeed();
   
   rnd.SetPrimesCouple(seed);
   double x = 0.; 

   double integral = 0.;
   int attempted = 0;
   int accepted = 0;

   double xNew = 0.;
   double A = 0.;
   
   double mu=1.;
   double sigma=1.;

   for (int i=0.; i<nstep; i++) {
      xNew = x + T(x);
      attempted ++;

      if(p(x,mu,sigma)==0.) A = 1.;
      else A = min(1.,p(xNew,mu,sigma)/p(x,mu,sigma));

      if (rnd.Rannyu() <= A){
         x = xNew;
         integral += x;
         accepted++;
      }
   }
   cout << "integrale di x: " <<  integral/accepted << endl;
   cout << "rate accettazione: " << double(accepted)/double(attempted) << endl;
   
   return 0;
}


double p(double x,double mu, double sigma) {   // funzione d'onda
   return exp(-pow((x-mu)/sigma,2)/2)+exp(-pow((x+mu)/sigma,2)/2);

}

double V(double x){
return pow(x,4)-5*pow(x,2)/2;
}

double T(double x) {
	return rnd.Rannyu() * delta - delta/2.;
}
