#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm> //serve per sort
#include <cmath>
#include "../random/random.h"
#include "../utils/utils.h"

using namespace std;

const int M_campionamenti = 100000;
const int N_blocchi = 100;
const int L_step = int(M_campionamenti/N_blocchi);
float dati[N_blocchi];
float medie[N_blocchi];
float deviazioni[N_blocchi];
Random rnd;
double delta = 1.6; 
double p(double x) { 
	if(x>0.) return exp(-2*x); 
	else return 0.; 
} 
double T(double x) { 
	return rnd.Rannyu() * delta - delta/2.; 
} 

double Integrand(double x) {
   return M_PI/2.*cos(M_PI/2.*x);
}
//Funzione peso per importance sampling
double Prob(double x) {
   return 2.*(1.-x);
}
void integraleMontecarlo()
{ 
  	vec r=zeros(M_campionamenti); 
  	vec s=zeros(M_campionamenti);
	double x = 0.; 
	double integral = 0.; 
	int attempted = 0; 
	int accepted = 0; 
	double xNew = 0.; 
	double A = 0.; 
	for (int i=0.; i<M_campionamenti; i++) { 
		xNew = x + T(x); 
		attempted ++; 
		if(p(x)==0.) A = 1.; 
		else A = min(1.,p(xNew)/p(x)); 
		if (rnd.Rannyu() <= A){ 
			x = xNew; 
			integral += x; 
			accepted++; 
		}
	       r(i) = Integrand(x); 
	       s(i) =  Prob(x) ;
	} 
	cout << "integrale di x: " <<  integral/accepted << endl; 
	cout << "rate accettazione: " 
		<< double(accepted)/double(attempted) << endl; 
   int N = 1e2;  // Number of blocks
   int L = M_campionamenti/N;  //# of numbers in a block
   mediaBlocchi2(r, N, L,"risultati/campionamentoUniforme.txt");
   mediaBlocchi2(s, N, L,"risultati/importanceSampling.txt");
}


int numCammini = 1e4;
int NBlocks = 1e2; // numero di blocchi in cui vengono raggruppati i random walks
int NSteps = 1e2; // passi cammino
void mediaBlocchiCammino(string);

//matrice di distanze al quadrato per ogni passo per ogni RW:
mat MSD(NSteps,  numCammini, arma::fill::zeros); 

void randomWalk(bool continuo) 
{
   // matrice che contiene la posizione raggiunta
   // da tutti i cammini a un certo istante:
   mat posizioniCasuali(3, numCammini, arma::fill::zeros);   
   
   std::ofstream filePosizioni;
   string nomeFilePos = "risultati/camminoDiscretoPos.txt";
   if(continuo) nomeFilePos = "risultati/camminoContinuoPos.txt";
   filePosizioni.open(nomeFilePos);
   filePosizioni.close();
   filePosizioni.open(nomeFilePos,ios::app);
   
   for (int l = 0; l < NSteps ; l++) 
   {
      for (int i = 0; i < numCammini; i++) 
      {
         if(!continuo)
         {
           int j = floor(rnd.Rannyu(0,3));
           if ( rnd.Rannyu() > 0.5 ) posizioniCasuali(j,i) += 1; 
           else posizioniCasuali(j,i) -= 1;
         } else
         {
           double theta = acos(1.-2.*rnd.Rannyu());
           double phi = 2*M_PI*rnd.Rannyu();
           posizioniCasuali(0,i) += sin(theta)*cos(phi);
           posizioniCasuali(1,i) += sin(theta)*sin(phi);
           posizioniCasuali(2,i) += cos(theta);
         }
         // salvo le posizioni finali raggiunte da ciascun cammino
         MSD(l,i) = pow(posizioniCasuali(0,i),2)+
         	pow(posizioniCasuali(1,i),2)+
         	pow(posizioniCasuali(2,i),2);
      }
      filePosizioni<<posizioniCasuali(0,0)<<"\t"
      <<posizioniCasuali(1,0)<<"\t"
      <<posizioniCasuali(2,0)<<"\t"<<endl;
   }
   filePosizioni.close();
   string risultati= "risultati/camminoDiscreto.txt"; 
   if(continuo) risultati= "risultati/camminoContinuo.txt"; 
   mediaBlocchiCammino(risultati);
}

void mediaBlocchiCammino(string risultati)
{
   vec SumAv(NBlocks,arma::fill::zeros);
   vec ErrAv(NBlocks,arma::fill::zeros);
   vec Sum2Av(NBlocks,arma::fill::zeros);
   mat Average(NSteps ,NBlocks); 
   //matrice 100x100. Le righe indicizzano il passo, 
   //le colonne la media sui blocchi del RW
   mat Average2(NSteps,NBlocks); 
   //numero di random walks per blocco
   int L = numCammini/NBlocks; 
   
   for (int l = 0; l < NSteps; l++) {
      for (int i = 0; i < NBlocks; i++) {
         double sum1 = 0.0;
         for (int j = 0; j < L; j++) {
            int k = j + i * L;
            sum1 += MSD(l,k);
         }
         //media sull'i-esimo blocco all'l-esimo passo
         Average(l,i) = sum1 / L; 
         Average2(l,i) = pow(Average(l,i),2);
      }
   }
   for (int l = 0; l < NSteps; l++) 
   {
        for (int i = 0; i <= NBlocks; i++) 
        {
            if(i<NBlocks){//////////////////////////////////////
            SumAv(l) += Average(l,i);
            Sum2Av(l) += Average2(l,i);}
        }
        SumAv(l)/=NBlocks;
        Sum2Av(l)/=NBlocks;
        ErrAv(l) = sqrt((Sum2Av(l) - SumAv(l)*SumAv(l))/NBlocks);
    }
    writeMeanError(sqrt(SumAv) ,  ErrAv/(2*sqrt(SumAv)) , risultati);                
}

int main ()
{
   rnd.SetSeed();   
   int seed=23; 
   rnd.SetPrimesCouple(seed);
   integraleMontecarlo();
   randomWalk(1);
   randomWalk(0);
   return 0;
}

