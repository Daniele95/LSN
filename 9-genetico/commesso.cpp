#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include "commesso.h"
#include <armadillo>

using namespace arma;
using namespace std;

float pswap=1;
float pshift=1;
float ppermut=0.5;
float pinversion=0.5;
float p = 5; //esponente per la selezione 
float pcrossover=0.5;

rowvec swap(rowvec);
bool check(mat);
double length(rowvec,int,mat);
double dist(colvec,colvec);
int Pbc(int j);
rowvec shift(rowvec path);

int populationSize = 500;
int Ncities = 20;
int generations = 600;
rowvec fittests(generations);

// array to store paths and their lengths
mat paths(populationSize, Ncities);
colvec lengths(populationSize);

// random cities positions
mat cities(2,Ncities); // 2 righe, N colonne

int main() {

	rnd.SetSeed();
   	rnd.SetPrimesCouple(23);
   	
	colvec city = linspace<colvec>(1,2,2);
	/*
	// generate random cities in square 1x1:
	for (int i=0;i<Ncities;i++) {
		city(0)=rnd.Rannyu();
		city(1)=rnd.Rannyu();
		cities.col(i)=city;
	}
	*/
	// generate random cities on circle:
	double angle;
	for (int i=0;i<Ncities;i++) {
		angle=rnd.Rannyu()*6.28;
		city(0)=cos(angle);
		city(1)=sin(angle);
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
		lengths(i) = length(newPath,Ncities,cities);
	}
	

	// sort paths with length
	uvec indices = sort_index(lengths);	
	paths = paths.rows(indices);
	lengths = lengths(indices);
	cout << "cammino più economico" << endl
	     << paths.row(0) << endl;

	cout << "evoluzione per " << generations
	     << " generazioni" << endl;

	for (int k = 0; k<generations; k++) {

		//mutazioni.
		//ciclo su tutti i cammini della generazione:
		for (int i=0; i<int(paths.n_rows);i++){
		  //muto con  probabilita p
		  if(rnd.Rannyu()<pswap)
		  	paths.row(i) = swap(paths.row(i));  
		  if(rnd.Rannyu()<pshift)
		  	paths.row(i) = shift(paths.row(i));
		  //ricalcolo la lunghezza:
		  lengths(i) = length(paths.row(i),
		    		      Ncities,cities);
		}
		
		// randomly select fit individuals
		// for reproduction
		// più alzo l'esponente, più vengono privilegiati 
		// i più fit
		int j = int(populationSize*pow(rnd.Rannyu(),p));
		
		// sort generation for fitness
		
		uvec indices = sort_index(lengths);	
		paths = paths.rows(indices);
		lengths = lengths(indices);	
		
		if( !check(paths) ) cout <<
			"qualche cammino non soddisfa"<<
			"le condizioni al contorno" << endl;
		
		fittests(k) = int(lengths(0));
	}
	colvec fitteststr=fittests.t();
	fitteststr.save("fittests.txt",raw_ascii);
	
	imat ipaths = conv_to<imat>::from(paths);
	ipaths.save("paths.txt",raw_ascii);
	
	ivec ilengths = conv_to<ivec>::from(lengths);
	ilengths.save("lengths.txt",raw_ascii);
	
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
e calcolo la funzione lunghezza con la distanza tra gli elementi
della matrice valutata alla posizione path[i]
*/

// costo di un cammino su un cerchio
double length(rowvec path,int Ncities, mat cities) {
	double length = 0.;
	for (int i = 1; i<Ncities; i++)
		length += dist(
		  cities.col(path[i]),cities.col(path[i-1])
		);
		
	length += dist(
	  cities.col(path[0]),cities.col(path[Ncities-1])
	);
	return length;
}
double dist(colvec city1, colvec city2){
	return pow((city1(0)-city2(0)),2)+
		pow((city1(1)-city2(1)),2);
}


int Pbc(int j) {
	if(j>Ncities-1) return j%Ncities+1;
	return j;
}
 
 
// Mutazioni-----------------------------------------------

/// 2) in posizione POS inverte M geni 
rowvec Inversion(int pos, int m,rowvec Gen){
   if(pos==0)       pos = 1;   
   if(pos>(Ncities-1))       pos = Pbc(pos);   
   if(m > Ncities-1)       m = 2;  

   int block[m];

   // metto da parte il blocco da swappare
   for(int j = 0; j < m; j++)      block[j] = Gen[Pbc(pos+j)];

   // riempio al contrario il buco lasciato dal blocco
   for(int i = 0; i<m; i++)      Gen[Pbc(pos+i)] = block[(m-1)-i];
	return Gen;
}

rowvec shift(rowvec path) {

	int m = floor(rnd.Rannyu(1,Ncities-1));
	for (int i=1; i<Ncities; i++) {
		rowvec pathTemp = path;
		path[i]=pathTemp[Pbc(i-m)];
	
	}
	return path;
}


/// 3) a partire da una posizione POS, sposta M geni adiacenti in avanti di N posizioni, 
///    (eccetto il primo gene e col vincolo m < Ncities-1)
rowvec Shift(int pos, int m, int n,rowvec Gen){
   if(pos==0)
      pos = 1;
   if(pos>(Ncities-1))   pos = Pbc(pos); 
   if(m >= Ncities-1)  m = 1;
   int block[m];
   // metto da parte il blocco da shiftare
   for(int j = 0; j < m; j++)  block[j] = Gen[Pbc(pos+j)];
   // shifto n volte i geni complementari per riempire il buco...
   for(int i = 0; i<n; i++)  Gen[Pbc(pos+i)] = Gen[Pbc(pos+m+i)];
   // ... e rimetto il blocco
   for(int i = 0; i<m; i++) Gen[Pbc(pos+n+i)] = block[i];
   return Gen;
}


/// 4) shifta tutto di "shift" posizioni
rowvec Shift2(int shift,rowvec Gen){
   
   int block[Ncities];

   // metto da parte il blocco da shiftare
   for(int j = 1; j < Ncities; j++)      block[j-1] = Gen[j];

   // sovrascrivo il blocco shiftato
   for(int i = 1; i < Ncities; i++)    Gen[Pbc(i+shift)] = block[i-1];
   return Gen;
}


/// 5) scambia M geni in posizione POS1 con altrettanti in POS2
rowvec MPermut(int pos1, int pos2, int m,rowvec Gen){
   if(pos1==0)      pos1 = 1;
   if(pos2==0)  pos2 = 1;
   if(pos1>(Ncities-1))  pos1 = Pbc(pos1);
   if(pos2>(Ncities-1)) pos2 = Pbc(pos2);
   if(m > Ncities/2)   m = 1;
   if(m > abs(pos2-pos1))  m = abs(pos1-pos2);
   if(pos1>pos2){
      int appo = pos2;
      pos2 = pos1;
      pos1 = appo;
   }

   int block[m];

   // metto da parte il blocco da swappare
   for(int j = 0; j < m; j++)   block[j] = Gen[Pbc(pos1+j)];

   // riempio il buco lasciato dal blocco...
   for(int i = 0; i<m; i++)  Gen[Pbc(pos1+i)] = Gen[Pbc(pos2+i)];
   
   // ... e rimetto il blocco in pos2
   for(int i = 0; i<m; i++)   Gen[Pbc(pos2+i)] = block[i];
   return Gen;
}


/*
/// 6) fa il crossover di len geni a partire dalla posizione pos, con un cromosoma parent2
void Crossover(int pos, int len, Chromosome parent2){

   if(pos==0) pos = 1;

   int block[len];
   // 1. cut their paths at the same position:
   //    metto da parte i blocchi da swappare
   //    (conservando la prima parte)
   for(int j = 0; j < len; j++) block[j] = Gen[Pbc(pos+j)];

   // 2. complete the paths with the missing cities adding them in the **order** 
   //    in which they appear in the consort (vale anche se sono separati!):
   int count1 = 0;
   
   for(int j = 1; j < Ncities; j++)
   
      // ciclo sui geni nei blocchi 
      for(int i = 0; i < len; i++)
         // quando trovo nell'altro genitore il gene presente
         // nel blocco, lo salvo nel buco scavato all'inizio
         if(parent2.Gen[Pbc(j)] == block[i] ){
            Gen[Pbc(pos+count1)] = parent2.Gen[Pbc(j)];
            count1++;
         }
      
}
*/
