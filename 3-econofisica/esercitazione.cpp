
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
      double Zi = rnd.Gauss(0,1);
      for (int j=0; j<nTempi; j++)
      {
         t[j+1] = (float(j)+1.)*T/float(nTempi);

         float z = rnd.Gauss(0.,1.);
         float dt = t[j+1]-t[j];
         S_[j+1]  = S[j] * exp((r-pow(sigma,2)/2.)*dt+sigma*z*sqrt(dt));
         
         // altre operazioni
         SpesaD(i) = SpesaD(i)*exp((Mean-1./2.*sigma*sigma)
         	*Incremento + sigma*Zi*Incremento);
         
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
     
   //cout << C << endl<<P<<endl;
     
     mediaBlocchi2( Co,  N, 
	 L,"risultati/outfileC.txt");
     mediaBlocchi2( Po,  N, 
	 L,"risultati/outfileP.txt");
     mediaBlocchi2( C2,  N, 
	 L,"risultati/outfileC2.txt");
     mediaBlocchi2( P2,  N, 
	 L,"risultati/outfileP2.txt");


}


double S(double S0, double m, double s, double t, double W) {
   return S0*exp((m-s*s/2.)*t+s*W);
}


int main ()
{
   rnd.SetSeed();   
   int seed=23; 
   rnd.SetPrimesCouple(seed);
   
   BlackScholes();


   return 0;
}

