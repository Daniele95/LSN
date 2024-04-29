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
#include <vector>
#include <cmath>
#include <algorithm> //serve per sort
#include "random.h"
#include "stat.h"

using namespace std;


//funzioni ricorrenti nell'esercizio.

//funzione integranda
double Integrand(double x) {
   return M_PI/2.*cos(M_PI/2.*x);
}

//
double Function(double x) {
   return M_PI*cos(M_PI*1./2.*x)/(4.*(1.-x));
}

//Funzione peso per importance sampling
double Prob(double x) {
   return 2.*(1.-x);
}

//Generatore di numeri distribuiti come la funzione Prob
double LinGen(double x) {
   return 1-sqrt(1-x);
}

void InitializeRand(Random& rnd); //Implementata sotto

 
int main (int argc, char *argv[]){

   Random rnd; /*dichiari la variabile rnd appartenente alla classe random, presa dai file random.cpp e random.h */
   InitializeRand(rnd);

   //ESERCIZIO 2.1.1

   int M = 1e6;              //Total number of throws
   int N = 1e2;                 // Number of blocks
   int L = M/N;    		//# of numbers in a block
   Stat st;

   //contiene i dati casuali grezzi
   vector<double> x(M);
   vector<double> r(M);
   for (int i = 0; i < M; i++) {
      x[i] = rnd.Rannyu();
      r[i] = Integrand(x[i]); // U[0,1) uniform distribution
   }

   //tengo salvati nel main solo i vettori contenenti i dati che poi verranno salvati su file. Gli altri creati e distrutti all'interno della funzione MeanAndErr	
   vector<double> sum_prog(N);
   vector<double> err_prog(N);
   fill(sum_prog.begin(), sum_prog.end(), 0.0);
   fill(err_prog.begin(), err_prog.end(), 0.0);

   st.MeanAndErr(r, sum_prog, err_prog, N, L); //media a blocchi e incertezze

    ofstream outfile211("output/outfile211.txt");
    for (int i = 0; i < N; ++i) {
       outfile211 << sum_prog[i] << "\t" << err_prog[i] << endl;
    }
    outfile211.close();


   //ESERCIZIO 2.1.2

   fill(sum_prog.begin(), sum_prog.end(), 0.0);
   fill(err_prog.begin(), err_prog.end(), 0.0);

   vector<double> s(M);
   for (int i = 0; i < M; i++) {
        s[i] = Function( LinGen(x[i]) ); // U[0,1) uniform distribution
   }

   st.MeanAndErr(s, sum_prog, err_prog, N, L);  //media a blocchi e incertezze

    ofstream outfile212("output/outfile212.txt");
    for (int i = 0; i < N; ++i) {
       outfile212 << sum_prog[i] << "\t" << err_prog[i] << endl;
    }
    outfile212.close();


   //ESERCIZIO 2.2.1

   int NRW = 1e4; //numero di random walks
   int NBlocks = 1e2; //numero di blocchi in cui vengono raggruppati i random walks
   L = NRW/NBlocks; //numero di random walks per blocco
   int NSteps = 1e2; //numero massimo di passi che si vogliono considerare per un RW


   vector<vector<int>> RW(3, vector<int>(NRW, 0));   //Matrice 3x1e4 posta a zero. Le righe inidicizzano la coordinata (x,y,z), le colonne il RW a passi fissati
   vector<vector<double>> Dist_Rad(NSteps, vector<double>(NRW, 0)); //matrice 1e2x1e4. Le righe indicizzano il passo, le colonne il RW

   //Calcolo le distanze al quadrato per ogni passo per ogni RW
   for (int l = 0; l < NSteps ; l++) { //ciclo su tutti i passi dei RW
      for (int i = 0; i < NRW; i++) { //ciclo su tutti i 1e4 RW
         int j = floor(rnd.Rannyu(0,3));
         if ( rnd.Rannyu() > 0.5 ) { RW[j][i] += 1; } else RW[j][i] -= 1;
         Dist_Rad[l][i] = pow(RW[0][i],2) + pow(RW[1][i],2) + pow(RW[2][i],2) ;          
      }
   }

   //Media a blocchi e incertezze
   vector<vector<double>> Average(NSteps, vector<double>(NBlocks)); //matrice 1e2x1e2. Le righe indicizzano il passo, le colonne la media sui blocchi del RW
   vector<vector<double>> Average2(NSteps, vector<double>(NBlocks));

   vector<double> SumAv(NBlocks);
   vector<double> ErrAv(NBlocks);
   vector<double> Sum2Av(NBlocks);
   fill(SumAv.begin(), SumAv.end(), 0.0);
   fill(ErrAv.begin(), ErrAv.end(), 0.0);

   for (int l = 0; l < NSteps; l++) {
      for (int i = 0; i < NBlocks; i++) {
         double sum1 = 0.0;
         for (int j = 0; j < L; j++) {
            int k = j + i * L;
            sum1 += Dist_Rad[l][k];
         }
         Average[l][i] = sum1 / L; //media sull'i-esimo blocco all'l-esimo passo
         Average2[l][i] = pow(Average[l][i], 2);
      }
   }

   for (int l = 0; l < NSteps; l++) {
        for (int i = 0; i <= NBlocks; i++) {
            SumAv[l] += Average[l][i];
            Sum2Av[l] += Average2[l][i];
        }
        SumAv[l]/=NBlocks;
        Sum2Av[l]/=NBlocks;
        ErrAv[l] = sqrt((Sum2Av[l] - SumAv[l]*SumAv[l])/NBlocks);
    }                     

    //Stampa
    ofstream outfile221("output/outfile221.txt");
    for (int l = 0; l < NSteps; ++l) {
       outfile221 << sqrt( SumAv[l] ) << "\t" << ErrAv[l]/(2*sqrt(SumAv[l])) << endl;
    }
    outfile221.close();


   //ESERCIZIO 2.2.2

   vector<vector<double>> RWC(3, vector<double>(NRW, 0));   //Versione continua
   vector<vector<double>> Dist_RadC(NSteps, vector<double>(NRW, 0)); //Versione continua

   //Calcolo le distanze al quadrato per ogni passo per ogni RW
   for (int l = 0; l < NSteps ; l++) { //ciclo su tutti i passi dei RW
      for (int i = 0; i < NRW; i++) { //ciclo su tutti i 1e4 RW
         double theta = acos(1.-2.*rnd.Rannyu());
         double phi = 2*M_PI*rnd.Rannyu();
         RWC[0][i] += sin(theta)*cos(phi);
         RWC[1][i] += sin(theta)*sin(phi);
         RWC[2][i] += cos(theta);
         Dist_RadC[l][i] = pow(RWC[0][i],2) + pow(RWC[1][i],2) + pow(RWC[2][i],2) ;          
      }
   }


   //Media a blocchi e incertezze
   fill(Sum2Av.begin(), Sum2Av.end(), 0.0);
   fill(SumAv.begin(), SumAv.end(), 0.0);
   fill(ErrAv.begin(), ErrAv.end(), 0.0);

   for (int l = 0; l < NSteps; l++) {
      for (int i = 0; i < NBlocks; i++) {
         double sum1 = 0.0;
         for (int j = 0; j < L; j++) {
            int k = j + i * L;
            sum1 += Dist_RadC[l][k];
         }
         Average[l][i] = sum1 / L; //media sull'i-esimo blocco all'l-esimo passo
         Average2[l][i] = pow(Average[l][i], 2);
      }
   }

   for (int l = 0; l < NSteps; l++) {
      for (int i = 0; i <= NBlocks; i++) {
         SumAv[l] += Average[l][i];
         Sum2Av[l] += Average2[l][i];
      }
      SumAv[l]/=NBlocks;
      Sum2Av[l]/=NBlocks;
      ErrAv[l] = sqrt((Sum2Av[l] - SumAv[l]*SumAv[l])/NBlocks);
    }                     

    //Stampa
    ofstream outfile222("output/outfile222.txt");
    for (int l = 0; l < NSteps; ++l) {
       outfile222 << sqrt( SumAv[l] ) << "\t" << ErrAv[l]/(2*sqrt(SumAv[l])) << endl;
    }
    outfile222.close();




   rnd.SaveSeed(); /*cambia i semi a fine lavoro*/

   return 0;
}




// ------------------------------//
// --------- FUNZIONI -----------//
// ------------------------------//

void InitializeRand(Random & rnd) {

   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ; //li legge dal file
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property; //variabile stringa chiamata property
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];//dai i semi all'algoritmo, presi dal file "seed.in"
            rnd.SetRandom(seed,p1,p2); 
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;
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
