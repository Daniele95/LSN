
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include "../random/random.h"
#include "../utils/utils.h"

using namespace std;
Random rnd;
void calcolaIntegrale() 
{

   int M = 1e6;              //Total number of throws
   int N = 1e2;                 // Number of blocks
   int L = M/N;    		//# of numbers in a block

   vec r(M);
   for (int i = 0; i < M; i++)  r(i) = rnd.Rannyu(); 
   mediaBlocchi2(r,N,L,"risultati/rMedia.txt");

   vec s(M);
   for(int i=0; i<M; i++)  s(i) = pow(r(i)-0.5,2);
   mediaBlocchi2(s, N, L,"risultati/rErrore.txt");


   //chi quadro

   int Nchi=1000; 
   int Lnew = M/Nchi;	
   mat matrix(N, Nchi,fill::zeros);
   for (int j = 0; j < Nchi; j++) 
   {
      vec VtoSort=r.subvec(j * Lnew, (j + 1) * Lnew - 1);
      VtoSort = sort(VtoSort);
      for (int i = 0; i < N; i++) {
         double count = 0;
         for (int k = 0; k < Lnew; k++) {
            if ((i / double(N) <= VtoSort(k)) && (VtoSort(k) < (i + 1) / double(N))) {
                count++;
            }
         }
         matrix(i,j) = count;
      }
   }

   vec ChiVect(Nchi);
   for (int j = 0; j < Nchi; j++) 
   {
      double sum = 0.;
      for (int i = 0; i < N; i++) 
          sum += pow(matrix(i,j) - double(Lnew/N), 2) / double(double(Lnew)/double(N));
      
      ChiVect(j) = sum;
   }

   ofstream outfileChi2("risultati/chiQuadro.txt");
   for (int i = 0; i < Nchi; ++i) outfileChi2 << ChiVect(i) << endl;
   outfileChi2.close();

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

