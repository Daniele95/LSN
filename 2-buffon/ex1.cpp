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
#include "../genRandom/random.h"

using namespace std;

const int M = 100000;
const int N = 100;
const int L = int(M/N);

float dati[N];
float medie[N];
float deviazioni[N];


float tiroCasuali() {

   Random rnd;
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

   ////////////////////////////////


   float d = 1.2; // spaziatura linee
   float l = 1.; // lunghezza ago

   float o=0.,x=0.,y=0.,r=0.,puntaAgo=0.;

   int tiri = 0;
   int intersezioni = 0;
   
   
   for(int i=0; i<N; i++)
   {
      
      while (tiri<L)
      {
         o = rnd.Rannyu()*d;
         // punto in un quadrato di raggio 2L
         x = (2.*rnd.Rannyu()-1.)*l;
         y = (2.*rnd.Rannyu()-1.)*l;
         
         r = sqrt(pow(x,2)+pow(y,2));
         
         puntaAgo = o+(x/r)*l;
         
         if( r<l )
         {
            tiri++;
            if (puntaAgo<0. || puntaAgo>d)
               dati[i] ++;               
         }
      
      }
      
      dati[i]  /= float(L);
      
      if(i==0) medie[i] = dati[i];
      
      for (int j=0; j<i; j++) {
         medie[i] += dati[j];
         if(i>0) deviazioni[i] += pow(dati[j],2);
      }
      if(i>0) {
         medie[i] /= float(i);
         deviazioni[i] = deviazioni[i]/float(i) - pow(medie[i],2);
         deviazioni[i] = sqrt(deviazioni[i] / float( i));
      }
      
      
      tiri = 0.;
   }
   
   cout << medie[89]<<endl;
   cout <<deviazioni[89]<<endl;
   float P = medie[99];
   rnd.SaveSeed();
   return (2.*l)/(P*d) ;
   
}


void write (float* array, int arrayLength, string nomeFile)
{
   std::ofstream ofile;
   ofile.open(nomeFile);

   for(int i=0; i<arrayLength; i++)
      ofile << array[i] << std::endl;

   ofile.close();   
}


int main (int argc, char *argv[]){


   cout << tiroCasuali() << endl;
   
   write(medie, N, "rMedia.txt");
   write(deviazioni, N, "rErrore.txt");
   
   
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
