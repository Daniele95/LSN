#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <armadillo>
#include "../random/random.h"


int seed[4];
Random rnd;

// accettazione metropolis
double centro=1.5; 
int accettazioni = 0;
Random& rnd;
int N_campionamenti=1e5;
int M_blocchi= 100;
int L_dimBlocco=N_campionamenti/M_blocchi;

rowvec posizioni(N_campionamenti); 
rowvec energia(N_campionamenti);
rowvec energiaMediaBlocchi(M_blocchi);
rowvec energiaMediaQuadraBlocchi(M_blocchi);
rowvec energiaErroreBlocchi(M_blocchi);

using namespace std;
using namespace arma;

// p(x) = |psi(x)|^2
double probability(double x, double m, double s) 
{
   double psi = exp(-0.5*pow((x-m)/s,2)+
   	 	exp(-0.5*pow((x+m)/s,2);
   return 0.282095*pow(psi,2)/s/(1+exp(-pow(m/s,2)));
}

// energia della particella alla pos. x
double energy(double x, double m, double s)
{
   double arg = m*x/(s*s);
   double epot=pow(x,4)-2.5*pow(x,2);
   double ekin = 
   	+0.5*(s*s-m*m)/pow(s,4)
   	-0.5*/pow(s,4)*pow(x,2)
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
        posizioni[i] = centro;
        energia[i] = energy(centro, media, sigma);
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
    x_proposto = rand.Gauss(x, passo);
    double alpha = 
    	probability(x_proposto, media, sigma) /
    	probability(x, media, sigma);
    if (rand.Rannyu() < alpha) {
        x = x_proposto;
        accettato = true;
    }
    return accettato;
}



// Simulated annealing algorithm using Metropolis block averages on a Gaussian probability distribution
void simulatedAnnealing(double mean, double stddev, int n_samples, int n_blocks, const string& energyFileName, const string& trajectoryFileName) {
    ofstream energyFile(energyFileName);
    ofstream trajectoryFile(trajectoryFileName);

    double temperature = 100.0; // Initial temperature
    double coolingRate = 0.003; // Cooling rate

    vec temperatures = linspace(temperature, 1.0, temperature / coolingRate + 1); 
    
   // Store the trajectory of mean and energy
    vec trajectory(n_samples); 

    double bestSolution = mean; // Initial solution
    double bestCost = valutaEnergia(mean, stddev, n_samples, n_blocks); // Initial cost

    energyFile << "# Temperature Steps Energy" << endl;

    for (double temp : temperatures) {
        double newMean = bestSolution + randu() * 0.1 - 0.05; // Generate a new solution near the current one
        double newCost = valutaEnergia(newMean, stddev, n_samples, n_blocks); // Calculate the cost of the new solution

        // If the new solution is better or satisfies the Boltzmann probability condition, accept it
        if (newCost < bestCost || exp(-(newCost - bestCost) / temp) > randu()) {
            bestSolution = newMean;
            bestCost = newCost;
        }

        // Write energy plot data to file
        energyFile << temp << " " << n_samples << " " << bestCost << endl;

        // Store trajectory data
        trajectory = join_horiz(trajectory, ones(n_samples) * temp);
        trajectory = join_horiz(trajectory, ones(n_samples) * bestCost);
    }

    energyFile.close();

    // Write trajectory data to file
    trajectoryFile << "# Temperature Energy" << endl;
    trajectory.save(trajectoryFile, raw_ascii);

    trajectoryFile.close();
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
   	sqrt((mediaQuadra(i) - media(i)*media(i))/n);
}

int main (int argc, char *argv[]){


   int seed=23;
   cout << "seme: "<<endl;
  
   rnd.SetSeed();
   
   
   mi[0]=15.5; sigma[0]=20.1;
   Calculate_H( sigma[0], mi[0], 1e4, 100); //1e4, 100 sono valori un po' bassi per il datablocking, velocizzano il codice

   string starting_mi = to_string(mi[0]);//variabile che definisce il nome del file di output

   H_SA[0] = H_sum_prog[n_blk-1];

   double beta_0 = beta;

   //ESERCIZIO 8.2.1 e 8.2.2--------------------------------------------------------
   for(int i=1; i<n_step; i++) { //ciclo sulle possibili temperature con cui fare SA

      beta += Delta_beta;
 
      delta_mi = rnd.Rannyu(-Lmi,Lmi)/pow(beta,1./4.); // legge di potenza con riduzione molto lenta. Anche se si parte da lontano l'esplorazione è maggiore
      delta_sigma = rnd.Rannyu(-Lsigma,Lsigma)/pow(beta,1./4.);
  
      mi[i] = mi[i-1] + delta_mi;
      sigma[i] = sigma[i-1] + delta_sigma;
      
      double H_old = H_sum_prog[n_blk-1]; //Variabile di supporto;
      double H_err_old = H_err_prog[n_blk-1]; //Variabile di supporto;

      Reset();
      Calculate_H(sigma[i],mi[i],n_samples,n_blk);

      double delta_H = H_sum_prog[n_blk-1] - H_old;
      double q = exp(-beta*(delta_H)); 
  
      cout << "mi    " << mi[i] << "       sigma  " << sigma[i] << endl;
 
      //Metropolis con SA
      if(rnd.Rannyu() < q) {
         H_SA[i] = H_sum_prog[n_blk-1];
         H_err_SA[i] = H_err_prog[n_blk-1];
      } else { //in caso di non accettazione si torna ai parametri iniziali
         H_SA[i] = H_old;
         H_err_SA[i] = H_err_old;
         mi[i] = mi[i]-delta_mi;
         sigma[i] = sigma[i]-delta_sigma;
      } 

   }

   ofstream outfile1("output_ESAeq.txt");
   const int wd=20;

   for (int i=0; i<n_step; i++) {
      outfile1 << beta_0+i*Delta_beta << setw(wd) << H_SA[i] << setw(wd) << H_err_SA[i] << setw(wd) << mi[i] << setw(wd) << sigma[i] << endl;
   }
   outfile1.close();

   //ESERCIZIO 8.2.3-----------------------------------------------------------
   //Una volta raggiunti valori stabili di mi e sigma (equilibrazione) si può campionare in modo più fine

   int step = 1500; //Deve essere <= di n_step!!!

   mi[0]=mi[n_step-1]; sigma[0]=sigma[n_step-1]; //i valori con cui si riparte sono quelli di equilibrazione
   Calculate_H( sigma[0], mi[0], 1e4, 100); //1e4, 100 sono valori un po' bassi per il datablocking, velocizzano il codice
   H_SA[0] = H_sum_prog[n_blk-1];

   beta_0=beta; //si riparte dalla temperatura di prima

   for(int i=1; i<step; i++) { //ciclo sulle possibili temperature con cui fare SA

      beta += Delta_beta;

      delta_mi = rnd.Rannyu(-Lmi,Lmi)/pow(beta,1./2.); // legge di potenza con riduzione più veloce
      delta_sigma = rnd.Rannyu(-Lsigma,Lsigma)/pow(beta,1./2.);
  
      mi[i] = mi[i-1] + delta_mi;
      sigma[i] = sigma[i-1] + delta_sigma;
      
      double H_old = H_sum_prog[n_blk-1]; //Variabile di supporto;
      double H_err_old = H_err_prog[n_blk-1]; //Variabile di supporto;

      Reset();
      Calculate_H(sigma[i],mi[i],n_samples,n_blk);

      double delta_H = H_sum_prog[n_blk-1] - H_old;
      double q = exp(-beta*(delta_H)); 
 
      cout << "mi    " << mi[i] << "       sigma  " << sigma[i] << endl;
 
      //Metropolis con SA
      if(rnd.Rannyu() < q) {
         H_SA[i] = H_sum_prog[n_blk-1];
         H_err_SA[i] = H_err_prog[n_blk-1];
      } else { //in caso di non accettazione si torna ai parametri iniziali
         H_SA[i] = H_old;
         H_err_SA[i] = H_err_old;
         mi[i] = mi[i]-delta_mi;
         sigma[i] = sigma[i]-delta_sigma;
      } 
 
   }
   ofstream outfile2("output_ESA.txt");
 
   for (int i=0; i<step; i++) {
      outfile2 << beta_0+i*Delta_beta << setw(wd) << H_SA[i] << setw(wd) << H_err_SA[i] << setw(wd) << mi[i] << setw(wd) << sigma[i] << endl;
   }
   outfile2.close();


   int Min_index = Minimum_Index(H_SA,step);
   ofstream outfile3("output_EnergyMin.txt");
   Reset();
   Calculate_H(sigma[Min_index],mi[Min_index],n_samples,n_blk);

   for (int i=0; i<n_blk; i++) {
      outfile3 << i << setw(wd) << H_sum_prog[i] << setw(wd) << H_err_prog[i] << setw(wd) << endl;
   }
   outfile3.close();


   //ESERCIZIO 8.2.4---------------------------------------------------------

   ofstream outfile4("output_Pos.txt");

   for (int i=0; i<n_samples; i++) {
      outfile4 << r[i] << endl;
   }
   outfile4.close();

  
   rnd.SaveSeed(); /*cambia i semi a fine lavoro*/

   return 0;
}


