
#include <iostream>
#include <fstream>
#include <cmath>
#include <functional> //serve per definire funzione di funzione
#include <iomanip> //serve per setprecision
#include "main.h"


using namespace std;



//Calcola HPsi/Psi
double H_Psi_over_PsiII(double x, double s, double m) {
    double DDPsi_over_Psi = ( m*m - 2*m*x*tanh(m*x/(s*s)) - s*s + x*x)/(2*pow(s,4));
    return -DDPsi_over_Psi + Potential(x);
}

//INTEGRALE DELL'ENERGIA-------------------------------------------------
//Calcola i valori dell'energia, la loro media a blocchi e il loro errore
void Calculate_H(double s,double m, int N_samp, int N_blocks) {

  // CICLO METROPOLIS---------------------------
   for (int i=0; i<50; i++) { //Ciclo affinchè i numeri siano generati con una distribuzione che si avvicina a quella asintotica
      MetrGauss(x, 1.8, rnd, s, m);
   }

   for (int i=0; i<N_samp; i++) {       
      if ( MetrGauss(x, 1.8, rnd, s, m) ) {Acceptance++;} //Modifica il salto in modo da avere accettazione intorno al 50%
      r[i]=x;
      H[i] = H_Psi_over_PsiII(x,s,m); //Calcolo delle energie
   }
   //cout << "Acceptance rate: " << double(Acceptance)/n_samples << endl;
   MeanAndErr(H,H_sum_prog,H_err_prog,N_blocks, N_samp/N_blocks );//media a blocchi e errore

}

//ELEMENTO DEL CICLO METROPOLIS
bool MetrGauss(double &x, double step_size, Random &rand, double s, double m) { 

    double x_proposed;
    bool A = false;

    x_proposed = rand.Gauss(x, step_size);

    double alpha = pow(Psi(x_proposed,s,m),2) / pow(Psi(x,s,m),2);

    if (rand.Rannyu() < alpha)
    {
        x = x_proposed;
        A = true;
    }
  
    return A;
}

int main (int argc, char *argv[]){

   Rnd_Inizialization();

   mi[0]=15.5; sigma[0]=20.1; //valori di partenza arbitrari
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



//----------------------------------------------------------------//

//----------------------------------------------------------------//

//VARIOUS FUNCTIONS THAT ALLOW TO COMPUTE THE HYDROGEN WAVEFUNCTION

//Calcola la funzione Psi
double Psi(double x, double s, double m) {
    double psi = exp( -pow( (x-m)/s , 2 )/2. ) + exp( -pow( (x+m)/s , 2 )/2. ) ;
    return psi;
}

//Calcola la derivata seconda di Psi
double DDPsi(double x, double s, double m) {
    double a1 = pow(x-m,2)*exp( -pow((x-m)/s,2) ) + pow(x+m,2)*exp( -pow((x+m)/s,2) );
    double DDPsi = 2./pow(s,2) * ( 2./pow(s,2) * a1 - Psi(x,s,m) );
    return DDPsi;
}

//Calcola il potenziale
double Potential(double x) {
    return pow(x,4)-5./2.*pow(x,2);
}



/*FONDAMENTALE IMPORTANTISSIMO: la variabile rnd deve essere passata per riferimento, altrimenti dopo che ha generato numero casuale non vengono aggiornati i parametri
interni e rigenera sempre lo stesso numero casuale!!
*/


void Reset() {

    Acceptance=0;
    for(int i=0; i<n_blk; i++) {
      H_sum_prog[i]=0;
      H_err_prog[i]=0;
    }
}

int Minimum_Index(double arr[], int size) {

    double minVal = arr[0]; // Assume the first element as the minimum value
    int index=0;
    for (int i = 1; i < size; i++) {
        if (arr[i] < minVal) {
            minVal = arr[i];
            index=i;
        }
    }

    return index;
}

double error(vector<double> &AV, vector<double> &AV2, int n) {
   if (n == 0) {
        return 0.0;
   } else {
        return sqrt((AV2[n] - AV[n]*AV[n])/n);
   }
}

void MeanAndErr(vector<double> &r, vector<double> &Mean, vector<double> &Errors, int N, int L) {

   vector<double> ave(N);
   vector<double> av2(N);
   vector<double> su2_prog(N);

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
        Mean[i]/=(i+1); //realizza il /n nel calcolo della media sui valori dei blocchi
        su2_prog[i]/=(i+1);
        Errors[i] = error(Mean, su2_prog, i);
    }
    return;
}
