#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include "commesso.h"
#include <armadillo>
#include <unordered_map>

using namespace arma;
using namespace std;


class cammino {

};
class popolazione {
    int N;

    void displayInfo() {
    }
};


float pswap=1;
float pshift=1;
float ppermut=0.5;
float pinversion=0.5;
float avversitaAmbientale = 5; //esponente per la selezione 
float pcrossover=0.5;

bool check(mat);
void ordinaCammini();
double lunghezzaCammino(rowvec cammino,mat mappa);
double dist(colvec,colvec);
int Pbc(int j);
void generateMap();

rowvec swap(rowvec);
rowvec trasponi_sottosequenza(rowvec);
rowvec scambia_sottosequenze(rowvec);
rowvec inverti_sottosequenza(rowvec);
mat riproduci(rowvec, rowvec);

int popolazione = 10;
int Ncities = 10;
int generations = 30;
rowvec fittests(generations);

// una generazione è la popolazione di cammini a un certo tempo,
// quindi una lista di cammini
mat generazione(popolazione, Ncities);
colvec lunghezze(popolazione);

// random cities positions
mat mappa(Ncities,2); // N righe,2 colonne, 

int main() {

	rnd.SetSeed();
   	rnd.SetPrimesCouple(23);
   	
   	generateMap();
	
	cout << "evolvo " << popolazione << " cammini "
	     << " di lunghezza " << Ncities << endl;

	// model path
	rowvec path = linspace<rowvec>(0,Ncities-1,Ncities);

	// popola la generazione iniziale di cammini
	for (int i=0;i<popolazione;i++) {
		rowvec newPath = swap(path);
		generazione.row(i) = newPath;
		lunghezze(i) = lunghezzaCammino(newPath,mappa);
	}
	
	generazione = generazione.rows(sort_index(lunghezze));
	//lunghezze = sort(lunghezze);
	
	cout << "cammino più economico" << endl
	     << generazione.row(0) << endl;

	cout << "evoluzione per " << generations
	     << " generazioni" << endl;

   // ciclo sulle generazioni
   for (int k = 0; k<generations; k++) {

      // creo la generazione successiva: seleziono
      // i cammini più brevi della generazione attuale
      mat nuovaGenerazione = generazione;
      
      
      // li muto 
      // e questi vanno a formare la generazione successiva
      
	// riempio la nuova generazione:
      for (int i=0; i<popolazione;i++){
         rowvec nuovoNato =
            generazione.row(int(popolazione*
               pow(rnd.Rannyu(),avversitaAmbientale)));
         //muto con  probabilita p
         if(rnd.Rannyu()<pswap) 
            generazione.row(i) = swap(generazione.row(i));  
         if(rnd.Rannyu()<pshift) 
            generazione.row(i) = inverti_sottosequenza(generazione.row(i));
		  
       
         nuovaGenerazione.row(i)=generazione.row(i);
      }
      generazione = nuovaGenerazione;
		
      // ordino la generazione in base alla lunghezza
      for (int i=0; i<popolazione;i++) 
         lunghezze(i) = lunghezzaCammino(generazione.row(i),mappa);
      generazione = generazione.rows(sort_index(lunghezze));
		
      if( !check(generazione) ) cout <<
      "qualche cammino non soddisfa"<<
      "le condizioni al contorno" << endl;
		
      fittests(k) = int(lunghezze(0));
   }
	colvec fitteststr=fittests.t();
	fitteststr.save("fittests.txt",raw_ascii);
	
	imat ipaths = conv_to<imat>::from(generazione);
	ipaths.save("paths.txt",raw_ascii);
	
	ivec ilengths = conv_to<ivec>::from(lunghezze);
	ilengths.save("lengths.txt",raw_ascii);
	
	cout << "cammino più economico" << endl
	     << generazione.row(0)<< endl;
	ivec iPath = conv_to<ivec>::from(generazione.row(0));
	iPath.save("Path.txt",raw_ascii);
	
	mat mappaOrdinata=mappa.rows(sort_index(generazione.row(0)));
	//mappaOrdinata.save("cities1.txt",raw_ascii);
	
	return 0;

}

bool check(mat generazione) {

	bool result = true;
	for (int i=0; i<popolazione; i++) { 
		rowvec riga = generazione.row(i);
		rowvec 	unici = unique(riga);
		result = ( unici.n_elem == riga.n_elem);
	}
	return result;

}

double lunghezzaCammino(rowvec cammino, mat mappa) {

   double lunghezza = 0.;
   for (int i = 1; i<Ncities; i++) // mappa.row(cammino[i]) è
   	// l'i-esima città visitata dal cammino
      lunghezza += norm(mappa.row(cammino[i])-mappa.row(cammino[i-1]));
   lunghezza += norm(mappa.row(cammino[0])-mappa.row(cammino[Ncities-1]));
   return lunghezza;
   
}

int Pbc(int j) {
	if(j>Ncities-1) return j%Ncities+1;
	return j;
}
 
int randInt(int a, int b){
   return floor(rnd.Rannyu(a,b));
}

void generateMap(){
	/*
	// generate random cities in square 1x1:
	for (int i=0;i<Ncities;i++) {
		citta(0)=rnd.Rannyu();
		cicittaty(1)=rnd.Rannyu();
		cities.col(i)=citta;
	}
	*/
	// generate random cities on circle:
	double angle;
	rowvec citta (2);
	for (int i=0;i<Ncities;i++) {
		angle=rnd.Rannyu()*6.28;
		citta(0)=cos(angle);
		citta(1)=sin(angle);
		mappa.row(i)=citta;
	}

	mappa.save("cities0.txt",raw_ascii);
	
}
//-------------------------------------------------------------
// MUTATIONS: they have to keep the first city in its place
// and the mutated sequence must still be a path

// swap two cities

rowvec swap(rowvec path) {

	path.swap_cols( randInt(1,Ncities), randInt(1,Ncities) );
	return path;
	
}

// taglia una sottosequenza 
// in corrispondenza di un certo locus
// e la incolla altrove

rowvec trasponi_sottosequenza(rowvec cammino) {

   int locus= randInt(1,Ncities-1);
   int lunghezza= randInt(1,Ncities-locus-1);
   int destinazione= randInt(1,Ncities-lunghezza-1);
   rowvec sottosequenza = 
   	cammino.subvec(locus, locus + lunghezza - 1);
   cammino.shed_cols(locus, locus + lunghezza - 1);
   cammino.insert_cols(destinazione, sottosequenza);
   return cammino;
   
}

rowvec scambia_sottosequenze(rowvec cammino) {

    int locus1 = randInt(1,Ncities-1);
    int lunghezza = randInt(1, (Ncities - locus1)/2);
    int locus2 = randInt(locus1+lunghezza,Ncities-lunghezza-1);
    rowvec sottosequenza1 = cammino.subvec(locus1, locus1 + lunghezza - 1);
    rowvec sottosequenza2 = cammino.subvec(locus2, locus2 + lunghezza - 1);
    cammino.subvec(locus1, locus1 + lunghezza - 1) = sottosequenza2;
    cammino.subvec(locus2, locus2 + lunghezza - 1) = sottosequenza1;
    return cammino;
    
}

rowvec inverti_sottosequenza(rowvec cammino) {

    int inizio= randInt(1,Ncities-1);
    int fine = randInt(inizio+1,Ncities);
    rowvec sottosequenza = cammino.subvec(inizio, fine - 1);
    reverse(sottosequenza.begin(), sottosequenza.end());
    cammino.subvec(inizio, fine - 1) = sottosequenza;
    return cammino;
    
}


mat riproduci(rowvec camminoPadre, rowvec camminoMadre) {

    int position = randInt(1, Ncities - 1);

    // Create a map to store the indices of numbers in camminoPadre
    unordered_map<double, int> indiciPadre;
    for (int i = 0; i < Ncities; ++i) indiciPadre[camminoPadre(i)] = i;

    // Rearrange the second part of camminoMadre according to the order in camminoPadre
    for (int i = position; i < Ncities; ++i) {
        auto it = indiciPadre.find(camminoMadre(i));
        if (it != indiciPadre.end()) {
            int index = it->second;
            if (index != i) {
                // Swap numbers
                double temp = camminoMadre(i);
                camminoMadre(i) = camminoMadre(index);
                camminoMadre(index) = temp;
            }
        }
    }

    // Rearrange the second part of camminoPadre according to the order in camminoMadre
    unordered_map<double, int> indiciMadre;
    indiciMadre.clear();
    for (int i = 0; i < Ncities; ++i) indiciMadre[camminoMadre(i)] = i;

    for (int i = position; i < Ncities; ++i) {
        auto it = indiciMadre.find(camminoPadre(i));
        if (it != indiciMadre.end()) {
            int index = it->second;
            if (index != i) {
                // Swap numbers
                double temp = camminoPadre(i);
                camminoPadre(i) = camminoPadre(index);
                camminoPadre(index) = temp;
            }
        }
    }
    return join_rows(camminoPadre, camminoMadre);
}
