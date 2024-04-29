
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


float tiroCasuali() {

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

   ////////////////////////////////


   float d = 1.2; // spaziatura linee
   float l = 1.; // lunghezza ago

   float o=0.,x=0.,y=0.,r=0.,puntaAgo=0.;

   int tiri = 0;
   int intersezioni = 0;
   
   
   for(int i=0; i<N_blocchi; i++)
   {
      
      while (tiri<L_step)
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
      
      dati[i]  /= float(L_step);
      
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
   //rnd.SaveSeed();
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





//------------- INTEGRAZIONE MONTE CARLO

int nstep=3000000; 
double delta = 1.6; 

double p(double x) { 
	if(x>0.) return exp(-2*x); 
	else return 0.; 
} 


double T(double x) { 
	return rnd.Rannyu() * delta - delta/2.; 
} 



void integraleMontecarlo(){ 

	double x = 0.; 

	double integral = 0.; 
	int attempted = 0; 
	int accepted = 0; 


	double xNew = 0.; 
	double A = 0.; 

	for (int i=0.; i<nstep; i++) { 
		xNew = x + T(x); 
		attempted ++; 

		if(p(x)==0.) A = 1.; 
		else A = min(1.,p(xNew)/p(x)); 

		if (rnd.Rannyu() <= A){ 
			x = xNew; 
			integral += x; 
			accepted++; 
		} 
	 
	} 
	cout << "integrale di x: " <<  integral/accepted << endl; 
	cout << "rate accettazione: " << double(accepted)/double(attempted) << endl; 

}



//-------------RANDOM WALK

   

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

 
void integraleMontecarlo2()
{

   int M = 1e6;              //Total number of throws
   int N = 1e2;                 // Number of blocks
   int L = M/N;    		//# of numbers in a block

   //contiene i dati casuali grezzi
   vec x(M);
   vec r(M);
   for (int i = 0; i < M; i++) {
      x(i) = rnd.Rannyu();
      r(i) = Integrand(x(i)); // U[0,1) uniform distribution
   }

   //tengo salvati nel main solo i vettori 
   //costd::ntenenti i dati che poi verranno 
   //salvati su file. Gli altri creati e 
   //distrutti all'interno della funzione MeanAndErr	
   vec sum_prog(N,arma::fill::zeros);
   vec err_prog(N,arma::fill::zeros);

   mediaBlocchi(r, sum_prog, err_prog, N, L); //media a blocchi e incertezze

    ofstream outfile211("risultati/outfile211.txt");
    for (int i = 0; i < N; ++i) 
       outfile211 << sum_prog(i) << "\t" << err_prog(i) << endl;
    outfile211.close();
    
   sum_prog.fill(0.);
   err_prog.fill(0.);
   
   vec s(M);
   for (int i = 0; i < M; i++) 
        s(i) = Function( LinGen(x(i)) ); 
        // U[0,1) uniform distribution
   
    mediaBlocchi(s, sum_prog, err_prog, N, L); 

    ofstream outfile212("risultati/outfile212.txt");
    for (int i = 0; i < N; ++i) 
       outfile212 << sum_prog(i) << "\t" << err_prog(i) << endl;
    outfile212.close();

    
} 




void randomWalkDiscreto()
{
   int NRW = 1e4; //numero di random walks
   int NBlocks = 1e2; //numero di blocchi in cui vengono raggruppati i random walks
   int L = NRW/NBlocks; //numero di random walks per blocco
   int NSteps = 1e2; //numero massimo di passi che si vogliono considerare per un RW
   
   vec SumAv(NBlocks,arma::fill::zeros);
   vec ErrAv(NBlocks,arma::fill::zeros);
   vec Sum2Av(NBlocks,arma::fill::zeros);

  //Media a blocchi e incertezze
   mat Average(NSteps ,NBlocks); 
   //matrice 1e2x1e2. Le righe indicizzano il passo, 
   //le colonne la media sui blocchi del RW
   mat Average2(NSteps,NBlocks); 

   mat RW(3, NRW, arma::fill::zeros);   
   //Matrice 3x1e4 posta a zero. Le righe
   // inidicizzano la coordinata (x,y,z), 
   //le colonne il RW a passi fissati
   mat Dist_Rad(NSteps,  NRW, arma::fill::zeros); 
   //matrice 1e2x1e4. Le righe indicizzano 
   //il passo, le colonne il RW
   //Calcolo le distanze al quadrato per ogni passo per ogni RW
   for (int l = 0; l < NSteps ; l++) { //ciclo su tutti i passi dei RW
      for (int i = 0; i < NRW; i++) { //ciclo su tutti i 1e4 RW
         int j = floor(rnd.Rannyu(0,3));
         if ( rnd.Rannyu() > 0.5 ) { RW(j,i) += 1; } else RW(j,i) -= 1;
         Dist_Rad(l,i) = pow(RW(0,i),2) + pow(RW(1,i),2) + pow(RW(2,i),2) ;          
      }
   }

   for (int l = 0; l < NSteps; l++) {
      for (int i = 0; i < NBlocks; i++) {
         double sum1 = 0.0;
         for (int j = 0; j < L; j++) {
            int k = j + i * L;
            sum1 += Dist_Rad(l,k);
         }
         Average(l,i) = sum1 / L; //media sull'i-esimo blocco all'l-esimo passo
         Average2(l,i) = pow(Average(l,i),2);
      }
   }
    
   for (int l = 0; l < NSteps; l++) {
        for (int i = 0; i <= NBlocks; i++) {
        
            if(i<NBlocks){//////////////////////////////////////
            SumAv(l) += Average(l,i);
            Sum2Av(l) += Average2(l,i);}
        }
        SumAv(l)/=NBlocks;
        Sum2Av(l)/=NBlocks;
        ErrAv(l) = sqrt((Sum2Av(l) - SumAv(l)*SumAv(l))/NBlocks);
    }                     

    //Stampa
    ofstream outfile221("risultati/outfile221.txt");
    for (int l = 0; l < NSteps; ++l) 
       outfile221 << sqrt( SumAv(l) ) << "\t" << ErrAv(l)/(2*sqrt(SumAv(l))) << endl;
    outfile221.close();

}
 
void randomWalkContinuo()
{
   int NRW = 1e4; //numero di random walks
   int NBlocks = 1e2; //numero di blocchi in cui vengono raggruppati i random walks
   int L = NRW/NBlocks; //numero di random walks per blocco
   int NSteps = 1e2; //numero massimo di passi che si vogliono considerare per un RW
   
   vec SumAv(NBlocks,arma::fill::zeros);
   vec ErrAv(NBlocks,arma::fill::zeros);
   vec Sum2Av(NBlocks,arma::fill::zeros);

  //Media a blocchi e incertezze
   mat Average(NSteps ,NBlocks); 
   //matrice 1e2x1e2. Le righe indicizzano il passo, 
   //le colonne la media sui blocchi del RW
   mat Average2(NSteps,NBlocks); 

   mat RWC(3, NRW, arma::fill::zeros);   //Versione continua
   mat Dist_RadC(NSteps, NRW, arma::fill::zeros); //Versione continua

   //Calcolo le distanze al quadrato per ogni passo per ogni RW
   for (int l = 0; l < NSteps ; l++) { //ciclo su tutti i passi dei RW
      for (int i = 0; i < NRW; i++) { //ciclo su tutti i 1e4 RW
         double theta = acos(1.-2.*rnd.Rannyu());
         double phi = 2*M_PI*rnd.Rannyu();
         RWC(0,i) += sin(theta)*cos(phi);
         RWC(1,i) += sin(theta)*sin(phi);
         RWC(2,i) += cos(theta);
         Dist_RadC(l,i) = pow(RWC(0,i),2) + 
         	pow(RWC(1,i),2) + pow(RWC(2,i),2) ;          
      }
   }

   //Media a blocchi e incertezze
   Sum2Av.fill(0.);
   SumAv.fill(0.);
   ErrAv.fill(0.);

   for (int l = 0; l < NSteps; l++) {
      for (int i = 0; i < NBlocks; i++) {
         double sum1 = 0.0;
         for (int j = 0; j < L; j++) {
            int k = j + i * L;
            sum1 += Dist_RadC(l,k);
         }
         Average(l,i) = sum1 / L; //media sull'i-esimo blocco all'l-esimo passo
         Average2(l,i) = pow(Average(l,i), 2);
      }
   }

   for (int l = 0; l < NSteps; l++) {
      for (int i = 0; i <= NBlocks; i++) {
            if(i<NBlocks){//////////////////////////////////////
         SumAv(l) += Average(l,i);
         Sum2Av(l) += Average2(l,i);}
      }
      SumAv(l)/=NBlocks;
      Sum2Av(l)/=NBlocks;
      ErrAv(l) = sqrt((Sum2Av(l) - SumAv(l)*SumAv(l))/NBlocks);
    }                     

    //Stampa
    ofstream outfile222("risultati/outfile222.txt");
    for (int l = 0; l < NSteps; ++l) 
       outfile222 << sqrt( SumAv(l) ) << "\t" << ErrAv(l)/(2*sqrt(SumAv(l))) << endl;
    outfile222.close();

}


int main (int argc, char *argv[]){


   rnd.SetSeed();   
   int seed=23; 
   rnd.SetPrimesCouple(seed);
   cout << tiroCasuali() << endl;
   
   write(medie, N_blocchi, "risultati/rMedia.txt");
   write(deviazioni, N_blocchi, "risultati/rErrore.txt");
   
   integraleMontecarlo2();
   randomWalkDiscreto();
   randomWalkContinuo();
   
   
   return 0;
}

