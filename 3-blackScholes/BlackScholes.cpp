/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <algorithm>
#include "../genRandom/random.h"

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
float sigma=0.25; // volatility


Random rnd;
   
void initRandom() {

   
   int seed[4];
   int p1, p2;
   ifstream Primes("../genRandom/Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("../genRandom/seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;
   
}


void write (float* array, int arrayLength, string nomeFile)
{
   std::ofstream ofile;
   ofile.open(nomeFile);

   for(int i=0; i<arrayLength; i++)
      ofile << array[i] << std::endl;

   ofile.close();   
}


void BlackScholes() {



   float S[M];   // prezzo al tempo T
   float Z[M]; // casuali gaussiani

   float C=0.;
   float P = 0.;

   for (int i=0; i<M; i++)
   {
   
      Z[i] = rnd.Gauss(0.,T);
      
     
      S[i] = S0*exp( (r-pow(sigma,2)/2.)*T+sigma*Z[i]*sqrt(T));
      
      
      // ricorsivo sui tempi //////////////////////////
      int nTempi = 100;
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
         
      }
      S[i] = S_[nTempi];
      ///////////////////////////////////
      
        
      C += exp(-r*T)*max(float(0.),float(S[i]-K));
      P+=exp(-r*T)*max(-float(S[i]-K),float(0.));
    
   }
   
   C /= M;
   P /= M;
     
   cout << C << endl<<P<<endl;
     
   rnd.SaveSeed();


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


int main (int argc, char *argv[]){
   initRandom();
   BlackScholes();

   return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
