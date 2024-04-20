#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <armadillo>
#include "../random/random.h"
#include <functional> 
#include <iomanip> 
#include <vector>

using namespace arma;
using namespace std;

const double pi=3.1415927;

int seed[4];
Random rnd;

// Simulation
int nstep=3000000;
double delta = 1.6;

double p(double,double,double);
double T(double);
double V(double);

//medie a blocchi
int n_samples=1e5; //# campionamenti di energia
int n_blk=100; //blocchi di energia
int L=n_samples/n_blk;
rowvec r(n_samples);//posizioni
rowvec H(n_samples);
rowvec H_sum_prog(n_blk);
rowvec H_err_prog(n_blk);

//Simulated Annealing
const int n_step=1500; //numero di temperature per il SA
double H_SA[n_step];
double H_err_SA[n_step];
double Lmi=1.5, Lsigma=1.5;
double delta_mi, delta_sigma;
double mi[n_step], sigma[n_step];

//Parametri
double x=1.5; //inizializzazione casuale
int Acceptance=0;
double temp=1., beta=1./temp, Delta_beta=2.;

//functions
double integrand(double,double,double);
void media(rowvec&,rowvec&,rowvec&,int,int);
double errore(rowvec&,rowvec&, int);

bool metropolis(double&,double,Random&,double,double);
void energia(double,double, int, int);

//int Minimum_Index(double[],int);
//double DDPsi(double,double,double);
//void Reset(void);

// v. anche lezione 1 e 4 per data blocking, e lez. 5 per calcolo sampling funzione con monte carlo

int main() {

   int seed=23;
   cout << "seme: "<<endl;
   
   rnd.SetSeed();
   
   rnd.SetPrimesCouple(seed);
   double x = 0.; 

   double integral = 0.;
   int attempted = 0;
   int accepted = 0;

   double xNew = 0.;
   double A = 0.;
   
   double mu=1.;
   double dev=1.;

   for (int i=0.; i<nstep; i++) {
      xNew = x + T(x);
      attempted ++;

      if(p(x,mu,dev)==0.) A = 1.;
      else A = min(1.,p(xNew,mu,dev)/p(x,mu,dev));

      if (rnd.Rannyu() <= A){
         x = xNew;
         integral += x;
         accepted++;
      }
   }
   cout << "integrale di x: " 
   	<<  integral/accepted << endl;
   cout << "rate accettazione: " 
   	<< double(accepted)/double(attempted) << endl;
   
   //valori di partenza arbitrari
   mi[0]=15.5; sigma[0]=20.1; 
   //1e4, 100 sono valori un po' bassi 
   //per il datablocking, velocizzano il codice
   energia( sigma[0], mi[0], 1e4, 100);
   //variabile che definisce il nome del file di output
   string starting_mi = to_string(mi[0]);
   H_SA[0] = H_sum_prog[n_blk-1];
   double beta_0 = beta;

   ofstream outfile1("output_ESAeq.txt");
   const int wd=20;

   for (int i=0; i<n_step; i++) 
      outfile1 << beta_0+i*Delta_beta 
      	<< setw(wd) << H_SA[i] 
      	<< setw(wd) << H_err_SA[i] 
      	<< setw(wd) << mi[i] 
      	<< setw(wd) << sigma[i] << endl;
  
   outfile1.close();

   return 0;
}

//elemento del ciclo d'integrazione
bool metropolis(double &x, double step_size, 
	Random &rand, double s, double m) { 
    double x_proposed;
    bool A = false;
    x_proposed = rand.Gauss(x, step_size);
    double alpha = pow(p(x_proposed,s,m),2) 
    	/ pow(p(x,s,m),2);
    if (rand.Rannyu() < alpha){
        x = x_proposed;
        A = true;
    }
    return A;
}

// HPsi/Psi
double integrand(double x, double s, double m) {
    return -( m*m - 2*m*x*tanh(m*x/(s*s)) - s*s + x*x)/(2*pow(s,4)) + V(x);
}

void energia(double s,double m, int N_samp, int N_blocks) {
   for (int i=0; i<50; i++)  metropolis(x, 1.8, rnd, s, m);
   for (int i=0; i<N_samp; i++) {       
      //Modifica il salto in modo
      // da avere accettazione intorno al 50%
      if ( metropolis(x, 1.8, rnd, s, m) ) Acceptance++;
      r[i]=x;
      H[i] = integrand(x,s,m);
   }
   media(H,H_sum_prog,H_err_prog,N_blocks, 
   	N_samp/N_blocks );
}

void media(rowvec& r, rowvec& Mean, rowvec& Errors, 
	int N, int L) {

   rowvec ave(N);
   rowvec av2(N);
   rowvec su2_prog(N);

   for (int i = 0; i < N; i++) {
      double sum1 = 0.0;
      for (int j = 0; j < L; j++) {
         int k = j + i * L;
         sum1 += r[k];
      }
      ave[i] = sum1 / L;
      av2[i] = pow(ave[i], 2);
   }
   for (int i = 0; i < N; i++) {
        for (int j = 0; j <= i; j++) {
            Mean[i] += ave[j];
            su2_prog[i] += av2[j];
        }
        //realizza il /n nel calcolo
        // della media sui valori dei blocchi
        Mean[i]/=(i+1); 
        su2_prog[i]/=(i+1);
        Errors[i] = errore(Mean, su2_prog, i);
    }
    return;
}

double errore(rowvec& AV, rowvec& AV2, int n) {
   if (n == 0) {
        return 0.0;
   } else {
        return sqrt((AV2[n] - AV[n]*AV[n])/n);
   }
}

double p(double x,double m, double s) {   // funzione d'onda
   return exp(-pow((x-m)/s,2)/2)+exp(-pow((x+m)/s,2)/2);

}

double V(double x){
return pow(x,4)-5*pow(x,2)/2;
}

double T(double x) {
	return rnd.Rannyu() * delta - delta/2.;
}
