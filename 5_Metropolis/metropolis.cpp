#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <armadillo>
#include <iomanip> 
#include "../random/random.h"
#include "../utils/utils.h"

using namespace std;
using namespace arma;

int seed[4];
Random rnd;
vec posizione = {-1.,2.5,2.}; //{-1.,2.5,2.}; 

double passo= 1.4; 
/*double passoGaussiano=0.8; //0.8
if(State) { 
   passoUniforme=2.5; 
   passoGaussiano=1.7;
}*/

int M_campionamenti = 1e6;
int N_blocchi=200;
int L_dimBlocco=M_campionamenti/N_blocchi;
	
int accettazione = 0; 
   
int M_equilibrazione = 400;

const double a0 = 1; // Bohr radius

//FUNZIONI
double Ground(double r, double theta, double phi){
    double psi = pow(M_E,-r)/sqrt(M_PI);
    return psi*psi;
}

double Excited(double r, double theta, double phi){
    double psi = r*pow(M_E,-r/2)*cos(theta)*sqrt(2/M_PI)/8.;
    return psi*psi;
}

double min (double a, double b) {
   if (a<b) {return a;}
   else {return b;}
}

vec OrtToSpher(double x, double y, double z) {
    vec SphericalCoord(3);
    SphericalCoord(0) = sqrt(x*x + y*y + z*z );
    SphericalCoord(1) = acos( z/sqrt(x*x + y*y + z*z ) );
    SphericalCoord(2) = (signbit(y) ? -1 : 1) * 
    	acos(x/sqrt(x*x + y*y)); 
    return SphericalCoord;
}
vec OrtToSpher(vec x) {
    vec SphericalCoord(3);
    SphericalCoord(0) = sqrt(x(0)*x(0) 
    	+ x(1)*x(1) +x(2)*x(2) );
    SphericalCoord(1) = acos( x(2)/sqrt(x(0)*x(0) 
    	+ x(1)*x(1) + x(2)*x(2) ) );
    SphericalCoord(2) = (signbit(x(1)) ? -1 : 1) * acos(x(0)/
    	sqrt(x(0)*x(0) + x(1)*x(1))); 
    return SphericalCoord;
}


bool Metropolis(vec &x, double passo, Random &rnd, 
	bool passoUniforme, bool groundState) { 
    vec x_Sph = OrtToSpher(x); 
    vec x_proposed(3);
    bool A = false;

    for (int j = 0; j < 3; j++) 
    {
       if(!passoUniforme) x_proposed(j) = 
          rnd.Gauss(x(j),passo); 
       else x_proposed(j) = 
          rnd.Rannyu(x(j)-passo,x(j)+passo);
    }

    vec x_proposed_Sph = OrtToSpher(x_proposed);
    double alpha;

    if (!groundState) {
    	 alpha = Ground(x_proposed_Sph(0),
    	 	x_proposed_Sph(1),x_proposed_Sph(2)) /
    	 	 Ground(x_Sph(0),x_Sph(1),x_Sph(2));
    	  } //se 0 (falso) considera ground
    else {
    	 alpha = Excited(x_proposed_Sph(0),
    	 	x_proposed_Sph(1),x_proposed_Sph(2)) /
    	 	Excited(x_Sph(0),x_Sph(1),x_Sph(2));
    }
    if (rnd.Rannyu() < alpha)
    {
        x = x_proposed;
        A = true;
    }
  
    return A;
}

void calculateWavefunction( bool groundState, bool passoUniforme)
{
   string stato;
   if(groundState) { stato = "Eccitato"; }
   else { stato = "Fondamentale"; }

   string tipoPasso;
   if(passoUniforme) { tipoPasso = "Uniforme"; }
   else { tipoPasso = "Gaussiano"; }

   //EQUILIBRAZIONE 
   ofstream outfileEq
      ("risultati/equilibrazione"+stato+tipoPasso+".txt");
      
   accettazione=0;
   for (int i=0; i<M_equilibrazione; i++)
   { 
      if ( Metropolis
      	(posizione, passo, rnd, passoUniforme, groundState) ) 
      	 	accettazione++;
      outfileEq << fixed << setprecision(5) 
      	<< Norm(posizione) << endl;
   }
   outfileEq.close();
   
   cout << "accettazione "+ tipoPasso+" "
   	<< double(accettazione)/M_equilibrazione << endl;

   //INIZIO SIMULAZIONE
   accettazione = 0; 

   vec posizioni(M_campionamenti); 

   for (int i=0; i<M_campionamenti; i++)
   {
      if ( Metropolis
      	(posizione, passo, rnd, passoUniforme, groundState) ) 
      	 accettazione++;

      if( i%100 == 0 ) 
      {
         ofstream outfilePos
            ("risultati/posizione"+stato+tipoPasso+".txt",ios::app);
            
         // STAMPO <r>
         outfilePos 
            << posizione(0) << setw(20) 
            << posizione(1) << setw(20) 
            << posizione(2) << endl;
         
         outfilePos.close();
      }
      posizioni(i) = Norm(posizione);
   }

   cout << "accettazione " + tipoPasso+" "  
   << double(accettazione)/M_campionamenti << endl;

   vec sum_prog(N_blocchi);
   vec err_prog(N_blocchi);

   mediaBlocchi(posizioni, sum_prog, 
   	err_prog, N_blocchi, L_dimBlocco);

   ofstream outfileR
   	("risultati/rMedio"+stato+tipoPasso+".txt");
   for (int i = 0; i < N_blocchi; ++i)    
       outfileR << sum_prog(i) 
       << setw(20) << err_prog(i) << endl;
   outfileR.close();
}

int main() 
{
   rnd.SetSeed();   
   int seed=23; 
   rnd.SetPrimesCouple(seed);	
   
   // passo uniforme:
   calculateWavefunction(0,0);
   calculateWavefunction(1,0);
   // passo gaussiano:
   calculateWavefunction(0,1);
   calculateWavefunction(1,1);
   
   return 0;
}


