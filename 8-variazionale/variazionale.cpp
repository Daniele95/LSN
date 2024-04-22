#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <armadillo>
#include <iomanip> 
#include "../random/random.h"

using namespace std;
using namespace arma;

int seed[4];
Random rnd;

const int passiAnnealing=1500; 
vec mu(passiAnnealing);
vec sigma(passiAnnealing);

bool metropolis(double&, double, double, double) ;
void mediaBlocchi();
double error(vec, vec, int i);

// accettazione metropolis
double centro=1.5; 
int accettazioni = 0;
int N_campionamenti=1e5;
int M_blocchi= 100;
int L_dimBlocco=N_campionamenti/M_blocchi;

vec posizioni(N_campionamenti); 
vec energia(N_campionamenti);
vec energiaMediaBlocchi(M_blocchi);
vec energiaErroreBlocchi(M_blocchi);

// p(x) = |psi(x)|^2
double probability(double x, double m, double s) 
{
   double psi = exp(-0.5*pow((x-m)/s,2))+
   	 	exp(-0.5*pow((x+m)/s,2));
   return 0.282095*pow(psi,2)/s/(1+exp(-pow(m/s,2)));
}

// energia della particella alla pos. x
double energy(double x, double m, double s)
{
   double arg = m*x/(s*s);
   double epot=pow(x,4)-2.5*pow(x,2);
   double ekin = 
   	 0.5*(s*s-m*m)/pow(s,4)
   	-0.5*pow(x,2)/pow(s,4)
   	+arg*tanh(arg)/(s*s);
   return ekin + epot;
}

// media dell'energia su stato gaussiano:
void valutaEnergia(double media, double sigma) 
{
    for (int i = 0; i < 50; i++) 
        metropolis(centro, 1.8, media, sigma);
             
    for (int i = 0; i < N_campionamenti; i++) {
        if (metropolis(centro, 1.8, media, sigma)) 
            accettazioni++;        
        posizioni(i) = centro;
        energia(i) = energy(centro, media, sigma);
    }
    mediaBlocchi();
}

// metropolis (esploro "annusando"
// la pdf con un cammino casuale)
bool metropolis(double& x, double passo, 
	 double media, double sigma) 
{
    double x_proposto;
    bool accettato = false;
    x_proposto = rnd.Gauss(x, passo);
    double alpha = 
    	probability(x_proposto, media, sigma) /
    	probability(x, media, sigma);
    if (rnd.Rannyu() < alpha) {
        x = x_proposto;
        accettato = true;
    }
    return accettato;
}

void mediaBlocchi() 
{
   vec energiaMediaTemp(M_blocchi);
   vec energiaErroreTemp(M_blocchi);

   for (int i = 0; i < M_blocchi; i++) 
   {
      energiaMediaTemp(i) =
         sum(energia.subvec(
            i * L_dimBlocco, 
            (i + 1) * L_dimBlocco - 1)) 
            / L_dimBlocco;
   }

   for (int i = 0; i < M_blocchi; i++) {
        for (int j = 0; j <= i; j++) 
        {
            energiaMediaBlocchi(i) 
            	+= energiaMediaTemp(j);
            energiaErroreTemp(i) 
            	+= pow(energiaMediaTemp(j), 2);
        } 
        // divido per il numero
        // di esperimenti fatti:
        energiaMediaBlocchi(i)/=(i+1);
        energiaErroreTemp(i)/=(i+1);
        energiaErroreBlocchi(i) =
           error(energiaMediaBlocchi,
            energiaErroreTemp, i);
    }
    return;
}

double error(vec media, vec mediaQuadra, int i) 
{
   if (i == 0) return 0.0;
   else return 
   	sqrt((mediaQuadra(i) - media(i)*media(i)) / i);
}

int Acceptance=0;
void Reset() 
{
    Acceptance=0;
    for(int i=0; i<M_blocchi; i++) {
      energiaMediaBlocchi[i]=0;
      energiaErroreBlocchi[i]=0;
    }
}


vec energiaAnnealing(passiAnnealing);
vec erroreAnnealing(passiAnnealing);
double Temperatura=1., beta=1./Temperatura, delta_beta=2.;
double Lmu=1.5, Lsigma=1.5;
double delta_mu, delta_sigma;
int main ()
{
   rnd.SetSeed();
      
   mu(0)=15.5; sigma(0)=20.1;
   string starting_mu = to_string(mu(0));
   
   valutaEnergia( mu(0), sigma(0) );
   energiaAnnealing(0) = energiaMediaBlocchi(M_blocchi-1);

   double beta_0 = beta;
   for(int i=1; i<passiAnnealing; i++) 
   {
      beta += delta_beta; 
      delta_mu = 
      	rnd.Rannyu(-Lmu,Lmu)/pow(beta,1./4.);
      delta_sigma = 
      	rnd.Rannyu(-Lsigma,Lsigma)/pow(beta,1./4.);
  
      mu(i) = mu(i-1) + delta_mu;
      sigma(i) = sigma(i-1) + delta_sigma;
      
      double energiaTemp = energiaMediaBlocchi(M_blocchi-1); 
      double erroreTemp = energiaErroreBlocchi(M_blocchi-1);

      Reset(); valutaEnergia(mu(i),sigma(i));

      double deltaEnergia = 
      	energiaMediaBlocchi(M_blocchi-1) - energiaTemp;
      double q = exp(-beta*(deltaEnergia)); 
  
      //Metropolis con simulated annealing
      if(rnd.Rannyu() < q) 
      {
         energiaAnnealing(i) = energiaMediaBlocchi(M_blocchi-1);
         erroreAnnealing(i) = energiaErroreBlocchi(M_blocchi-1);
      } else 
      {
         energiaAnnealing(i) = energiaTemp;
         erroreAnnealing(i) = erroreTemp;
         mu(i) = mu(i) - delta_mu;
         sigma(i) = sigma(i) - delta_sigma;
      } 
   }

   ofstream outfile1(
   	"risultati/energiaAnnealingEquilibrio.txt");
   const int wd=20;

   for (int i=0; i < passiAnnealing; i++) {
      outfile1 << beta_0 + i * delta_beta << setw(wd) 
      	<< energiaAnnealing(i) << setw(wd) 
      	<< erroreAnnealing(i) << setw(wd) 
      	<< mu(i) << setw(wd) 
      	<< sigma(i) << endl;
   }
   outfile1.close();

   int passiAnnealing2 = passiAnnealing; 
   // valori di equilibrio:
   mu(0)=mu(passiAnnealing-1); 
   sigma(0)=sigma(passiAnnealing-1);
   
   valutaEnergia(  mu(0), sigma(0)); 
   energiaAnnealing(0) = energiaMediaBlocchi(M_blocchi-1);

   beta_0=beta; 
   
   for(int i=1; i < passiAnnealing2; i++) 
   {
      beta += delta_beta;
      delta_mu = 
         rnd.Rannyu(-Lmu,Lmu) / pow(beta,.5);       
      delta_sigma = 
      	 rnd.Rannyu(-Lsigma,Lsigma) / pow(beta,.5);
  
      mu(i) = mu(i-1) + delta_mu;
      sigma(i) = sigma(i-1) + delta_sigma;
      
      double energiaTemp = energiaMediaBlocchi(M_blocchi-1);
      double erroreTemp = energiaErroreBlocchi(M_blocchi-1);

      Reset(); valutaEnergia(mu(0), sigma(0));

      double deltaEnergia = 
      	energiaMediaBlocchi(M_blocchi-1) - energiaTemp;
      double q = exp(-beta*(deltaEnergia)); 
  
      //Metropolis con simulated annealing
      if(rnd.Rannyu() < q) 
      {
         energiaAnnealing(i) = energiaMediaBlocchi(M_blocchi-1);
         erroreAnnealing(i) = energiaErroreBlocchi(M_blocchi-1);
      } else 
      {
         energiaAnnealing(i) = energiaTemp;
         erroreAnnealing(i) = erroreTemp;
         mu(i) -= delta_mu;
         sigma(i) -= delta_sigma;
      } 
   }
   ofstream outfile2("risultati/energiaAnnealing.txt");
 
   for (int i=0; i<passiAnnealing2; i++) 
   {
      outfile2 << beta_0 + i * delta_beta << setw(wd) 
      << energiaAnnealing(i) << setw(wd) 
      << erroreAnnealing(i) << setw(wd) 
      << mu(i) << setw(wd) 
      << sigma(i) << endl;
   }
   outfile2.close();


   ofstream outfile3("risultati/energiaMinima.txt");   
   int indiceMin = energiaAnnealing.index_min();
   Reset();valutaEnergia(mu(indiceMin), sigma(indiceMin));   
   for (int i=0; i<M_blocchi; i++) 
   {
      outfile3 << i << setw(wd) 
      << energiaMediaBlocchi(i) << setw(wd) 
      << energiaErroreBlocchi(i) << setw(wd) 
      << endl;
   }
   outfile3.close();


   ofstream outfile4("risultati/istogramma-posizioni.txt");
   for (int i=0; i<N_campionamenti; i++)
       outfile4 << posizioni(i) << endl;   
   outfile4.close();  
   rnd.SaveSeed(); 
   return 0;
}


