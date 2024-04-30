
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include "../random/random.h"
#include "../utils/utils.h"

using namespace std;

const int M = 100000;
static float r[M];         

Random rnd;

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

void lorenziane() {

   int Nrep = 1e5;
   int dimension[] = {1, 2, 10, 100};
   double Mean = 1.;

   for(int l=0; l<4; l++) {

      vector<double> S_N(Nrep);
      for (int i=0; i<Nrep; i++) {
         for (int j=0; j<dimension[l]; j++) {
            S_N[i]+=rnd.Exp(Mean);
         }
         S_N[i]=S_N[i]/dimension[l];
      }

      string Index = to_string(dimension[l]);
      ofstream outfileLCTexp("risultati/outfileLCTexp"+Index+".txt");
      for (int i = 0; i < Nrep; ++i) {
           outfileLCTexp << S_N[i] << endl;
      }
      outfileLCTexp.close();
   }


   double G = 1.;
   for(int l=0; l<4; l++) {

      vector<double> S_N(Nrep);
      for (int i=0; i<Nrep; i++) {
         for (int j=0; j<dimension[l]; j++) {
            S_N[i]+=rnd.Lor(G);
         }
         S_N[i]=S_N[i]/dimension[l];
      }

      string Index = to_string(dimension[l]);
      ofstream outfileLCTlor("risultati/outfileLCTlor"+Index+".txt");
      for (int i = 0; i < Nrep; ++i) {
           outfileLCTlor << S_N[i] << endl;
      }
      outfileLCTlor.close();
   }

   double Length = 1.;
   double d = 2.;
   int M_campionamenti = 1e5;         //Total number of throws
   int N_blocchi = 1e2;         // Number of blocks
   int L_step = M_campionamenti/N_blocchi; 

   vector<double> b(M_campionamenti);
   for (int i = 0; i < M_campionamenti; i++) 
        b[i] = 1.*d*rnd.Rannyu(); 
       // U[0,1) uniform distribution
   

   vector<double> l(M_campionamenti);

   for (int i = 0; i < M_campionamenti; i++) 
      l[i] = Length*sin(2.*rnd.UnPhiAR());
   

   vector<double> ave_(N_blocchi);
   vector<double> av2_(N_blocchi);
   vector<double> sum_prog_(N_blocchi); //nuovi vettori
   vector<double> su2_prog_(N_blocchi);
   vector<double> err_prog_(N_blocchi);

    for (int i = 0; i < N_blocchi; i++) {
        int Nhit = 0;
        for (int j = 0; j < L_step; j++) {
            int k = j + i * L_step;
            if ( b[k]+l[k] < 0. || b[k]+l[k] > d ) {
                Nhit++;
            }
        }
        ave_[i] = 2.0 * Length * L_step / (d * Nhit); 
        // calculation of pi for the i-th block
        av2_[i] = pow(ave_[i], 2); // (r_i)^2
    }

   for (int i = 0; i < N_blocchi; i++) {
        for (int j = 0; j <= i; j++) {
            sum_prog_[i] += ave_[j];
            su2_prog_[i] += av2_[j];
        }
        sum_prog_[i] /=(i+1);
        su2_prog_[i] /=(i+1);
        err_prog_[i] = error(sum_prog_, su2_prog_, i);
    }

    ofstream outfile13("risultati/outfile13.txt");
    for (int i = 0; i < N_blocchi; ++i) {
        outfile13 << sum_prog_[i] << "\t" << err_prog_[i] << endl;

    }
    outfile13.close();
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


int main (int argc, char *argv[])
{
   int esercizio = 1;
   tiroCasuali();            // tiro M casuali

   rnd.SetSeed();   
   int seed=23; 
   rnd.SetPrimesCouple(seed);

   esperimenti( esercizio );         // riempio gli array  
   analisiDati();
  
   write(sum_prog, N, "risultati/rMedia.txt");
   write(err_prog, N, "risultati/rErrore.txt");
   
   
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
   
   lorenziane();
   
   return 0;
}

