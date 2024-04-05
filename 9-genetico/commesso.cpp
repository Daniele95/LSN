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
double cost(rowvec,int,mat);
double dist(colvec,colvec);
int Pbc(int j);
rowvec shift(rowvec path);

int populationSize = 10;
int Ncities = 15;
int generations = 1000;
rowvec fittests(generations);

// array to store paths and their cost
mat paths(populationSize, Ncities);
colvec costs(populationSize);

// random cities positions
mat cities(2,Ncities); // 2 righe, N colonne

int main() {

	rnd.SetSeed();
   	rnd.SetPrimesCouple(23);
	
	// generate random cities in square 1x1:
	colvec city = linspace<colvec>(1,2,2);
	for (int i=0;i<Ncities;i++) {
		city(0)=rnd.Rannyu();
		city(1)=rnd.Rannyu();
		cities.col(i)=city;
	}
	mat citiestr = cities.t(); //transpose, for read in python
	citiestr.save("cities0.txt",raw_ascii);
	
	cout << "evolvo " << populationSize << " cammini "
	     << " di lunghezza " << Ncities << endl;

	// model path
	rowvec path = linspace<rowvec>(0,Ncities-1,Ncities);

	// randomly create
	// population of paths
	for (int i=0;i<populationSize;i++) {
		rowvec newPath = swap(path);
		paths.row(i) = newPath;
		costs(i) = cost(newPath,Ncities,cities);
	}
	
	float p = 5;

	// sort paths with cost
	uvec indices = sort_index(costs);	
	paths = paths.rows(indices);
	costs = costs(indices);
	cout << "cammino più economico" << endl
	     << paths.row(0) << endl;

	cout << "evoluzione per " << generations
	     << " generazioni" << endl;

	for (int k = 0; k<generations; k++) {

		// select randomly the fittests	
		// più alzo l'esponente, più vengono privilegiati 
		// i più fit
		int j = int(populationSize*pow(rnd.Rannyu(),p));

		paths.row(j) = swap(paths.row(j));
		costs(j) = cost(paths.row(j),Ncities,cities);
		
		// sort
		uvec indices = sort_index(costs);	
		paths = paths.rows(indices);
		costs = costs(indices);	

		if( !check(paths) ) cout <<
			"qualche cammino non soddisfa"<<
			"le condizioni al contorno" << endl;
		
		fittests(k) = int(costs(0));
	}
	imat ipaths = conv_to<imat>::from(paths);
	ipaths.save("paths.txt",raw_ascii);
	
	ivec icosts = conv_to<ivec>::from(costs);
	icosts.save("costs.txt",raw_ascii);
	
	cout << "cammino più economico" << endl
	     << paths.row(0)<< endl;
	ivec iPath = conv_to<ivec>::from(paths.row(0));
	iPath.save("Path.txt",raw_ascii);
	
	uvec sorted = sort_index(paths.row(0));	
	mat citiessorted = cities.cols(sorted);	
	mat citiessortedtr = citiessorted.t();
	citiessortedtr.save("cities1.txt",raw_ascii);
	
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
	int index1 = floor(rnd.Rannyu(1,Ncities));
	int index2 = floor(rnd.Rannyu(1,Ncities));
	int temp = path[index1];
	path[index1] = path[index2];
	path[index2] = temp;
	return path;
}

/* creo matrice di città (x e y), 
la riempio con numeri casuali tra 0 e 1, 
e calcolo la funzione costo con la distanza tra gli elementi
della matrice valutata alla posizione path[i]
*/

// costo di un cammino su un cerchio
double cost(rowvec path,int Ncities, mat cities) {
	double cost = 0.;
	for (int i = 1; i<Ncities; i++)
		cost += dist(
		  cities.col(path[i]),cities.col(path[i-1])
		);
		
	cost += dist(
	  cities.col(path[0]),cities.col(path[Ncities-1])
	);
	return cost;
}
double dist(colvec city1, colvec city2){
	return pow((city1(0)-city2(0)),2)+
		pow((city1(1)-city2(1)),2);
}


rowvec shift(rowvec path) {

	int m = floor(rnd.Rannyu(1,Ncities-1));
	for (int i=1; i<Ncities; i++) {
		rowvec pathTemp = path;
		path[i]=pathTemp[Pbc(i-m)];
	
	}
	return path;
}


int Pbc(int j) {
	if(j>Ncities-1) return j%Ncities+1;
	return j;
}
 
