#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <armadillo>
#include "../random/random.h"

// accettazione metropolis
double centro=1.5; 
int accettazioni = 0;
Random& rnd;
int M_campionamenti=1e5;
int N_blocchi= 100;

rowvec posizioni(M_campionamenti); 
rowvec energia(M_campionamenti);
rowvec media_prog(N_blocchi);
rowvec errore_prog(N_blocchi);

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
             
    for (int i = 0; i < M_campionamenti; i++) {
        if (metropolis(centro, 1.8, media, sigma)) 
            accettazioni++;        
        posizioni[i] = centro;
        energia[i] = energy(centro, media, sigma);
    }

    mediaBlocco();
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

// Calcola la media a blocchi e l'errore
void mediaBlocco() {
    int M_campionamenti = dati.n_elem;
    int lunghezza_blocco = M_campionamenti / N_blocchi;
    vec media(N_blocchi), media_quad(N_blocchi), somma_prog(N_blocchi), somma_prog_quad(N_blocchi);

    for (int i = 0; i < N_blocchi; i++) {
        int inizio = i * lunghezza_blocco;
        int fine = (i + 1) * lunghezza_blocco;

        media[i] = mean(dati.subvec(inizio, fine - 1));
        media_quad[i] = pow(media[i], 2);

        if (i == 0) {
            somma_prog[i] = media[i];
            somma_prog_quad[i] = media_quad[i];
        } else {
            somma_prog[i] = somma_prog[i - 1] + media[i];
            somma_prog_quad[i] = somma_prog_quad[i - 1] + media_quad[i];
        }
    }

    for (int i = 0; i < N_blocchi; i++) {
        media_blocchi[i] = somma_prog[i] / (i + 1);
        errore_blocchi[i] = sqrt((somma_prog_quad[i] / (i + 1) - pow(media_blocchi[i], 2)) / i);
    }
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

int main() {
    double mean = 0.0; // Mean of the Gaussian distribution
    double stddev = 1.0; // Standard deviation of the Gaussian distribution
    int n_samples = 10000; // Number of Metropolis samples
    int n_blocks = 100; // Number of blocks for block averaging
    string energyFileName = "energy_plot.dat"; // File to write energy plot data
    string trajectoryFileName = "trajectory_plot.dat"; // File to write trajectory data

    simulatedAnnealing(mean, stddev, n_samples, n_blocks, energyFileName, trajectoryFileName);

    return 0;
}

