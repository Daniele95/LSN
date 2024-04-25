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

ofstream outFilePos;

int seed[4];
Random rnd;
vec posIniziale ={-1.,20.,2.};
vec posizione(3);

int M_campionamenti = 1e6;
int N_blocchi=200;
int L_dimBlocco=M_campionamenti/N_blocchi;
	
int accettazione = 0; 
   
int M_equilibrazione = 400;

double probab0(vec v)
{
	return pow(M_E,-2.*Norm(v))/M_PI;
}

double probab1(vec v)
{
	double r=Norm(v);
	double cosTheta=v(2)/r;
	double psi=r*pow(M_E,-r*0.5)*cosTheta;
	return psi*psi/(32*M_PI);
}

double min (double a, double b) {
   if (a<b) {return a;}
   else {return b;}
}

bool metropolis(vec &x, double passo, Random &rnd, 
	bool passoUniforme, bool groundState) 
{ 
   vec x_proposed={0.,0.,0.};
   bool A = false;
   if(!passoUniforme) x_proposed = rnd.Gauss(x,passo); 
   else x_proposed = rnd.Rannyu(x, passo);
   double alpha;

   if (!groundState) alpha = probab0(x_proposed)/probab0(x);
   else alpha = probab1(x_proposed)/probab1(x);
   
   if (rnd.Rannyu() < alpha)
   {
       x = x_proposed;
       A = true;
   }
   return A;
}

void scriviSuFile(vec3 posizione,string filePosizioni)
{
	outFilePos.open(filePosizioni,ios::app);   	
         outFilePos 
            << posizione(0) << setw(20) 
            << posizione(1) << setw(20) 
            << posizione(2) << endl;
         outFilePos.close();
}

void calculateWavefunction( bool groundState, 
	bool passoUniforme,double passo)
{
   string stato;
   if(groundState) stato = "Eccitato";
   else stato = "Fondamentale";

   string tipoPasso;
   if(passoUniforme) tipoPasso = "Uniforme";
   else tipoPasso = "Gaussiano"; 

   string filePosizioni=
      "risultati/posizione"+stato+tipoPasso+".txt";
   
   outFilePos.open(filePosizioni,ios::trunc);
   outFilePos.close();
   
   //EQUILIBRAZIONE -------------------------------------
   posizione=posIniziale;
   accettazione=0;
   for (int i=0; i<M_equilibrazione; i++)
   { 
      if ( metropolis(posizione, passo,
       rnd, passoUniforme, groundState) ) 
      	 	accettazione++;    
       scriviSuFile(posizione,filePosizioni);
   }
   
   cout << "accettazione "+ tipoPasso+" "
   	<< double(accettazione)/M_equilibrazione << endl;

   //INIZIO SIMULAZIONE -------------------------------------
   accettazione = 0; 

   vec posizioni(M_campionamenti); 

   for (int i=0; i<M_campionamenti; i++)
   {
      if ( metropolis
      	(posizione, passo, rnd, passoUniforme, groundState) ) 
      	 accettazione++;
      posizioni(i) = Norm(posizione);
       // STAMPO <r>
      if( i%100 == 0 ) scriviSuFile(posizione,filePosizioni);
   }

   cout << "accettazione " + tipoPasso+" "  
   << double(accettazione)/M_campionamenti << endl;

   vec sum_prog(N_blocchi);
   vec err_prog(N_blocchi);
   
   mediaBlocchi(posizioni, sum_prog, 
   	err_prog, N_blocchi, L_dimBlocco);
   
   string fileMedie=
   	"risultati/rMedio"+stato+tipoPasso+".txt";
   ofstream outFileBlockAverage(fileMedie);
   for (int i = 0; i < N_blocchi; ++i)    
       outFileBlockAverage << sum_prog(i) 
       << setw(20) << err_prog(i) << endl;
   outFileBlockAverage.close();
}

int main() 
{
   rnd.SetSeed();   
   int seed=23; 
   rnd.SetPrimesCouple(seed);	
   
   // passo uniforme:
   calculateWavefunction(0,0,1.4);
   calculateWavefunction(1,0,2.5);
   
   // passo gaussiano:
   calculateWavefunction(0,1,0.8);
   calculateWavefunction(1,1,1.7);
   return 0;
}


