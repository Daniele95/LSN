
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <algorithm>
#include "../random/random.h"
#include "../utils/utils.h"

using namespace std;

const int M = 100000;
const int N = 100;
const int L = int(M/N);

float dati[N];
float medie[N];
float deviazioni[N];


float S0=100.; // asset price at S(0)
float T=1.; // delivery time
float K=100.; // strike price
float r=0.1; // risk-free interest rate
float Mean=r; // risk-free interest rate
float sigma=0.25; // volatility


Random rnd;
   
double So(double S0, double m, double s, double t, double W) {
   return S0*exp((m-s*s/2.)*t+s*W);
}
double Max(double a, double b) {
   if (a>b) {return a;}
   else return b;
}
void BlackScholes() 
{
   float S[M];   // prezzo al tempo T
   float Z[M]; // casuali gaussiani

   float C=0.;
   float P = 0.;

// altri vettori
   vec W(M);
   vec Spesa(M);
   vec Co(M);
   vec Po(M);
   vec C2(M);
   vec P2(M);
   vec SpesaD(M);

   for (int i=0; i<M; i++)
   {
      SpesaD[i] = S0;
      //
      Z[i] = rnd.Gauss(0.,T);
      S[i] = S0*exp( (r-pow(sigma,2)/2.)*T+sigma*Z[i]*sqrt(T));
      // ricorsivo sui tempi //////////////////////////
      int nTempi = 100;
      double Incremento = T/nTempi;
      float S_[nTempi];
      S_[0] = S0;
      float t[nTempi];
      t[0] = 0.;
      for (int j=0; j<nTempi; j++)
      {
         t[j+1] = (float(j)+1.)*T/float(nTempi);

         float z = rnd.Gauss(0.,1.);
         float dt = t[j+1]-t[j];
         S_[j+1]  = S[j] * exp((r-pow(sigma,2)/2.)*dt+sigma*z*sqrt(dt));
         
         // altre operazioni
         SpesaD(i) = SpesaD(i)*exp((Mean-1./2.*sigma*sigma)
         	*Incremento + sigma*z*Incremento);
         
      }
      S[i] = S_[nTempi];
      ///////////////////////////////////
      
        
      C += exp(-r*T)*max(float(0.),float(S[i]-K));
      P+=exp(-r*T)*max(-float(S[i]-K),float(0.));
    
     // altre operazioni
     
      W(i) = rnd.Gauss(0,T);
      Spesa(i) = So(S0,Mean,sigma,T,W(i));
      Co(i) = exp(-Mean*T)*Max(0.,Spesa(i)-K);
      Po(i) = exp(-Mean*T)*Max(0.,K-Spesa(i));
      C2(i) = exp(-Mean*T)*Max(0.,SpesaD(i)-K);
      P2(i) = exp(-Mean*T)*Max(0.,K-SpesaD(i));
   }
   
   C /= M;
   P /= M;
     
   cout << C << endl<<P<<endl;
     mediaBlocchi2( Co,  N, 
	 L,"risultati/outfileC.txt");
     mediaBlocchi2( Po,  N, 
	 L,"risultati/outfileP.txt");
     mediaBlocchi2( C2,  N, 
	 L,"risultati/outfileC2.txt");
     mediaBlocchi2( P2,  N, 
	 L,"risultati/outfileP2.txt");
	
   //rnd.SaveSeed();


/* perché mi viene 
18.7554
0.0895789
 invece che i risultati giusti
  == BLACK-SCHOLES ==
call:  14.975790778311286
put:  5.4595325819072364
 ???
 */

}


double S(double S0, double m, double s, double t, double W) {
   return S0*exp((m-s*s/2.)*t+s*W);
}
/*

void blackScholes2() {


   //ESERCIZIO 3.1.1

   int M = 1e5;              //Total number of throws
   int NBlocks = 1e2;                 // Number of blocks
   int L = M/NBlocks;    		//# of numbers in a block

   //Parametri caratteristici
   double Mean = 0.1;
   double Sigma = 0.25;
   double S0 = 100.;
   double T = 1.;
   double k = 100.;


   vec W(M);
   for (int i=0; i<M; i++) {
      W(i) = rnd.Gauss(0,T);
   }   

   vec Spesa(M);
   for (int i=0; i<M; i++) {
      Spesa(i) = S(S0,Mean,Sigma,T,W(i));
   }   

   vec C(M);
   vec P(M);
   for (int i=0; i<M; i++) {
      C(i) = exp(-Mean*T)*Max(0.,Spesa(i)-k);
      P(i) = exp(-Mean*T)*Max(0.,k-Spesa(i));
   }  

   vec MeanC(NBlocks);
   vec MeanP(NBlocks);
   vec ErrorC(NBlocks);
   vec ErrorP(NBlocks);

   mediaBlocchi(C,MeanC,ErrorC,NBlocks,L);
   mediaBlocchi(P,MeanP,ErrorP,NBlocks,L);

   ofstream outfile311("risultati/outfile311.txt");
    for (int i = 0; i < NBlocks; ++i) 
       outfile311 << MeanC(i) << "\t" << ErrorC(i) 
       	<< "\t" << MeanP(i) << "\t" << ErrorP(i) << endl;
    outfile311.close();


   //ESERCIZIO 3.1.2

   int NIntervalli = 100;
   double Incremento = T/NIntervalli;
   
   vec SpesaD(M);
   for (int i=0; i<M; i++) {
      SpesaD[i] = S0;
   }      

   for (int i=0; i<M; i++) {
      double Zi = rnd.Gauss(0,1);
      for (int j=0; j<NIntervalli; j++) {
         SpesaD(i) = SpesaD(i)*exp((Mean-1./2.*Sigma*Sigma)
         	*Incremento + Sigma*Zi*Incremento);
      }   
   }   

   for (int i=0; i<M; i++) {
      C(i) = exp(-Mean*T)*Max(0.,SpesaD(i)-k);
      P(i) = exp(-Mean*T)*Max(0.,k-SpesaD(i));
   }  
   MeanC.fill(0.);
   MeanP.fill(0.);
   ErrorC.fill(0.);
   ErrorP.fill(0.);

   mediaBlocchi(C,MeanC,ErrorC,NBlocks,L);
   mediaBlocchi(P,MeanP,ErrorP,NBlocks,L);

   ofstream outfile312("risultati/outfile312.txt");
    for (int i = 0; i < NBlocks; ++i) 
       outfile312 << MeanC(i) << "\t" << ErrorC(i) 
       	<< "\t" << MeanP(i) << "\t" << ErrorP(i) << endl;
    outfile312.close();


}
*/
int main (int argc, char *argv[]){


   rnd.SetSeed();   
   int seed=23; 
   rnd.SetPrimesCouple(seed);
   
   BlackScholes();

  // blackScholes2();

   return 0;
}

