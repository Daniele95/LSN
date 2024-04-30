
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include "../random/random.h"
#include "../utils/utils.h"

using namespace std;
Random rnd;

const int M = 100000;
 float r[M];         

float integranda(float r_k, int esercizio)
{
   if ( esercizio == 1 ) return r_k;
   if ( esercizio == 2 ) return pow(r_k-0.5,2);
   else return 0.;
}

const int Niniziale=100;
const int Liniziale= int(M/Niniziale);
float ave[Niniziale];                // ciascun esperimento ha come output
float av2[Niniziale];                // una media e una varianza
float sum_prog[Niniziale];           // analisi dati: ogni nuovo esperimento
float su2_prog[Niniziale];           // calcolo media e deviazione standard
float err_prog[Niniziale];           // su tutti gli esperimenti passati
	
void esperimenti( int N, int L, int esercizio )
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

void analisiDati()
{
   for ( int i = 0; i < Niniziale; i ++ )
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

// N = 100, L = 1000
float chiSquare[Niniziale];
void chiSquareTest()
{
   for ( int i = 0; i < Niniziale; i++ )
   {
      int howManyInRange = 0;
      for ( int j = 0; j < Liniziale; j++ )
      {
         int k = j + i * Liniziale;
         // conto quanti tra i 1000 numeri casuali nel blocco i
         // ricadono tra i/100 e (i+1)/100: dovrebbero essere
         // 1000/100 = 10
         if ( ( r[ k ] > float(i) / Niniziale ) && ( r[ k ] < float( i + 1 ) / Niniziale ) )
            howManyInRange ++;
      }
      chiSquare[i] = pow( howManyInRange - Liniziale/Niniziale, 2) 
      	/ (Liniziale/Niniziale);
   }
}

void calcolaIntegrale()
{
   int N = Niniziale;
   int L = Liniziale;  // suddivido i M numeri casuali
   // in N esperimenti contenenti
   // ciascuno L numeri casuali
   int esercizio = 1;
   for(int i=0; i<M; i++) r[i] = rnd.Rannyu();
   esperimenti(N,L, esercizio );         // riempio gli array  
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
   esperimenti(N,L, esercizio );         
   analisiDati();
   write(sum_prog, N, "risultati/sigmaMedia.txt");
   write(err_prog, N, "risultati/sigmaErrore.txt");   
   chiSquareTest();
   write(chiSquare, N, "risultati/chiQuadro.txt");
   int sumChiSquare = 0;
   for ( int j = 0; j < N; j++ )
      sumChiSquare += chiSquare[j];      
   //cout <<"chi quadro = "<< sumChiSquare << endl;
}

void lorenziane() 
{
   int Nrep = 1e5;
   ivec dimension = {1, 2, 10, 100};
   double Mean = 1.;
   for(int l=0; l<4; l++) 
   {
      vec S_N(Nrep);
      for (int i=0; i<Nrep; i++) { 
         for (int j=0; j<dimension(l); j++) {
            S_N(i)+=rnd.Exp(Mean);
         }
         S_N(i)=S_N(i)/dimension(l);
      }
      string Index = to_string(dimension[l]);
      ofstream outfileLCTexp("risultati/outfileLCTexp"+Index+".txt");
      for (int i = 0; i < Nrep; ++i) {
           outfileLCTexp << S_N(i) << endl;
      }
      outfileLCTexp.close();
   }
   double G = 1.;
   for(int l=0; l<4; l++) 
   {
      vec S_N(Nrep);
      for (int i=0; i<Nrep; i++) {
         for (int j=0; j<dimension(l); j++) {
            S_N(i)+=rnd.Lor(G);
         }
         S_N(i)=S_N(i)/dimension(l);
      }
      string Index = to_string(dimension[l]);
      ofstream outfileLCTlor("risultati/outfileLCTlor"+Index+".txt");
      for (int i = 0; i < Nrep; ++i) {
           outfileLCTlor << S_N(i) << endl;
      }
      outfileLCTlor.close();
   }
}

float buffon() {
   int M_campionamenti = 100000;
   int N_blocchi = 100;
   int L_stepblocco = int(M_campionamenti/N_blocchi);
   vec dati(N_blocchi);
   vec medie(N_blocchi);
   vec deviazioni(N_blocchi);
   float d = 1.2; // spaziatura linee
   float l = 1.; // lunghezza ago
   float o=0.,x=0.,y=0.,r=0.,puntaAgo=0.;
   int tiri = 0;
   for(int i=0; i<N_blocchi; i++)
   {
      while (tiri<L_stepblocco)
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
               dati(i) ++;               
         }
      }
      dati(i)  /= float(L_stepblocco);
      if(i==0) medie[i] = dati(i);
      for (int j=0; j<i; j++) {
         medie(i) += dati(j);
         if(i>0) deviazioni(i) += pow(dati(j),2);
      }
      if(i>0) {
         medie[i] /= float(i);
         deviazioni(i) = deviazioni(i)/float(i) - pow(medie(i),2);
         deviazioni(i) = sqrt(deviazioni(i) / float( i));
      }
      tiri = 0.;
   }
   cout << medie[89]<<endl;
   cout <<deviazioni[89]<<endl;
   float P = medie[99];
   writeMeanError(medie, deviazioni, "risultati/buffon.txt");
   return (2.*l)/(P*d) ;
}

int main (int argc, char *argv[])
{
   rnd.SetSeed();   
   int seed=23; 
   rnd.SetPrimesCouple(seed);
   calcolaIntegrale();
   lorenziane();
   cout << buffon() << endl;
   return 0;
}

