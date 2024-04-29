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
#include "../random/random.h"

using namespace std;

const int M = 100000;
static float r[M];         

void tiroCasuali() {

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("../random/Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("../random/seed.in");
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

                                    // genero M numeri 
                                    // con una distribuzione uniforme tra 0 e 1
   for(int i=0; i<M; i++)
      r[i] = rnd.Rannyu();

   //rnd.SaveSeed();
   
}

const int N = 100;
const int L = int(M/N);             // suddivido i M numeri casuali
                                    // in N esperimenti contenenti
                                    // ciascuno L numeri casuali
static float ave[N];                // ciascun esperimento ha come output
static float av2[N];                // una media e una varianza

float integranda(float r_k, int esercizio)
{
   if ( esercizio == 1 ) return r_k;
   if ( esercizio == 2 ) return pow(r_k-0.5,2);
   else return 0.;
}


void esperimenti( int esercizio )
{
   
   for (int i = 0; i < N; i ++)     // in ciascun esperimento
   {
      float sum = 0;
      
      for (int j = 0; j < L; j ++)  // calcolo l'integrale
      {                             // su un gruppo di L numeri casuali
         int k = j + i * L;
         sum += integranda( r[k], esercizio );
      }
      
      ave[i] = sum / L;
      av2[i] = pow((ave[i]),2);
   }
   
}


float error( float* AV, float* AV2, int n )
{
   if ( n == 0 ) return 0.;
   else return sqrt(  ( AV2[ n ] - pow( AV[ n ], 2) ) / n );
}


static float sum_prog[N];           // analisi dati: ogni nuovo esperimento
static float su2_prog[N];           // calcolo media e deviazione standard
static float err_prog[N];           // su tutti gli esperimenti passati
                                    
void analisiDati()
{
   for ( int i = 0; i < N; i ++ )
   {
      for ( int j = 0; j < i + 1; j++ )
      {
         sum_prog[i] += ave[j];
         su2_prog[i] += av2[j];
      }
      sum_prog[i] /= i + 1;
      su2_prog[i] /= i + 1;
      err_prog[i] = error( sum_prog, su2_prog, i );
   }
}

void write (float* array, int arrayLength, string nomeFile)
{
   std::ofstream ofile;
   ofile.open(nomeFile);

   for(int i=0; i<arrayLength; i++)
      ofile << array[i] << std::endl;

   ofile.close();   
}

   
// N = 100, L = 1000
static float chiSquare[N];

void chiSquareTest()
{
   for ( int i = 0; i < N; i++ )
   {
      int howManyInRange = 0;
      for ( int j = 0; j < L; j++ )
      {
         int k = j + i * L;
         // conto quanti tra i 1000 numeri casuali nel blocco i
         // ricadono tra i/100 e (i+1)/100: dovrebbero essere
         // 1000/100 = 10
         if ( ( r[ k ] > float(i) / N ) && ( r[ k ] < float( i + 1 ) / N ) )
            howManyInRange ++;
      }
      chiSquare[i] = pow( howManyInRange - L/N, 2) / (L/N);
      
   }
}


int main (int argc, char *argv[]){

   int esercizio = 1;

   tiroCasuali();            // tiro M casuali

   esperimenti( esercizio );         // riempio gli array  
   analisiDati();
  
   write(sum_prog, N, "rMedia.txt");
   write(err_prog, N, "rErrore.txt");
   
   
   for ( int i = 0; i < N; i++ )     // resetto gli accumulatori
   {
      ave[i] = 0.;
      av2[i] = 0.;
      sum_prog[i] = 0.;
      su2_prog[i] = 0.;
      err_prog[i] = 0.;
   }
   
   esercizio = 2;
   esperimenti( esercizio );         
   analisiDati();
  
   write(sum_prog, N, "risultati/sigmaMedia.txt");
   write(err_prog, N, "risultati/sigmaErrore.txt");   
   
   
   chiSquareTest();
   write(chiSquare, N, "risultati/chiQuadro.txt");
   int sumChiSquare = 0;
   for ( int j = 0; j < N; j++ )
      sumChiSquare += chiSquare[j];      
   cout << sumChiSquare << endl;
   
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
