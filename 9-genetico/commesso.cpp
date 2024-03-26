#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include "commesso.h"
#include <armadillo>

using namespace arma;
using namespace std;

rowvec swap(rowvec);
bool check(mat);
double cost(rowvec);
int Pbc(int j);
rowvec shift(rowvec path);

int populationSize = 10;
int pathLength = 20;
int generations = 100;
rowvec fittests(generations);

int main() {

        // Read seed for random numbers
        int p1, p2;
        ifstream Primes("Primes");
        Primes >> p1 >> p2;
        Primes.close();

        ifstream input("seed.in");
        input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
        rnd.SetRandom( seed, p1, p2 );
        input.close();
	//

	cout << "evolvo " << populationSize << " cammini "
		<< " di lunghezza " << pathLength << endl;

	// model path
	rowvec path = linspace<rowvec>(1,pathLength,pathLength);
	// array to store paths and their cost
	mat paths(populationSize, pathLength);
	colvec costs(populationSize);

	// randomly create
	// population of paths
	for (int i=0;i<populationSize;i++) {
		rowvec newPath = swap(path);
		paths.row(i) = newPath;
		costs(i) = cost(newPath);
	}
	
	float p = 5;

	// sort paths with cost
	uvec indices = sort_index(costs);	
	paths = paths.rows(indices);
	costs = costs(indices);
	cout << paths.row(0) << endl;

	cout << "evoluzione per " << generations
	       	<< " generazioni" <<endl;

	for (int k = 0; k<generations; k++) {

		// select randomly the fittests	
		// più alzo l'esponente, più vengono privilegiati 
		// i più fit
		int j = int(populationSize*pow(rnd.Rannyu(),p));

		paths.row(j) = swap(paths.row(j));
		costs(j) = cost(paths.row(j));
		
		// sort
		uvec indices = sort_index(costs);	
		paths = paths.rows(indices);
		costs = costs(indices);	

		if( !check(paths) ) cout <<
			"qualche cammino non soddisfa"<<
			"le condizioni al contorno" << endl;
		fittests(k) = int(costs(0));
	}

	fittests.save("costs.txt",raw_ascii);

	cout << paths.row(0)<< endl;
	return 0;

}

bool check(mat paths) {

	bool result = true;
	for (int i=0; i<populationSize; i++) { 
		rowvec riga = paths.row(i);
		rowvec 	unici = unique(riga);
		result = ( unici.n_elem == riga.n_elem);
	}
	return result;

}


// swap two cities in the path
// keeping the first in place
rowvec swap(rowvec path) {
	int index1 = floor(rnd.Rannyu(1,pathLength));
	int index2 = floor(rnd.Rannyu(1,pathLength));
	int temp = path[index1];
	path[index1] = path[index2];
	path[index2] = temp;
	return path;
}


double cost(rowvec path) {
	double cost = 0.;
	for (int i = 1; i<pathLength; i++)
		cost += abs(path[i]-path[i-1]);
	cost += abs(path[0]-path[pathLength-1]);
	return cost;
}


rowvec shift(rowvec path) {
	

	int m = floor(rnd.Rannyu(1,pathLength-1));
	for (int i=1; i<pathLength; i++) {
		rowvec pathTemp = path;
		path[i]=pathTemp[Pbc(i-m)];
	
	}
	
}


int Pbc(int j) {
	if(j>4) return j%5+1;
	return j;
}
 
