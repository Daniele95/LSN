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


float avversitaAmbientale =3; //esponente per la selezione 
float pMutazione=1.;
float pScambio=0.1*pMutazione;
float pTrasponi=0.1*pMutazione;
float pScambiaSequenze=0.1*pMutazione;
float pInverti=0.1*pMutazione;
float pRiproduzione=0.5;

bool check(mat);
void ordinaCammini();
double lunghezzaCammino(rowvec cammino,mat mappa);
double dist(colvec,colvec);
int Pbc(int j);
void generaMappaCerchio();
void generaMappaQuadrato();

rowvec muta(rowvec);
rowvec swap(rowvec);
rowvec trasponi_sottosequenza(rowvec);
rowvec scambia_sottosequenze(rowvec);
rowvec inverti_sottosequenza(rowvec);
mat riproduci(rowvec, rowvec);

int popolazione = 300;
int Ncities_max=100;
int Ncities = 34;
int generazioni = 1000;
int generazioni_max = 2000;
colvec migliori(generazioni_max);
colvec migliori_semimedia(generazioni_max);
void risolviCommessoViaggiatore();

// una generazione è la popolazione di cammini a un certo tempo,
// quindi una lista di cammini.
// un cammino è un vettore riga. 
// per fare un array di cammini, li incolonno

mat generazione(popolazione, Ncities_max);
colvec lunghezze(popolazione);

// random cities positions
mat mappa(Ncities_max,2); // N righe,2 colonne, 
mat mappaOrdinata;

int main() {

   int seed;
   cout << "seme: "<<endl;
   cin >> seed;

   rnd.SetSeed();
   rnd.SetPrimesCouple(seed);
   	
   
   cin>> avversitaAmbientale; //esponente per la selezione 
   cin>> pScambio;
   cin>> pTrasponi;
   cin>> pScambiaSequenze;
   cin>> pInverti;
   cin>> pRiproduzione;
	
   cin>> popolazione;
   cin>> Ncities;
   cin>> generazioni;

   // taglio le matrici al valore letto per Ncities e generazioni
   mappa=mappa.rows(0,Ncities); // N righe,2 colonne, 
   generazione=generazione.cols(0,Ncities-1);	
   migliori=migliori.rows(0,generazioni-1);
   migliori_semimedia=migliori_semimedia.rows(0,generazioni-1);
   	
   generaMappaCerchio();
   risolviCommessoViaggiatore();
   mappaOrdinata.save("risultati/cerchio1.txt",raw_ascii);
   migliori.save("risultati/miglioriCerchio.txt",raw_ascii);   
   migliori_semimedia.save("risultati/migliori_semimediaCerchio.txt",raw_ascii);
   
   generaMappaQuadrato();
   risolviCommessoViaggiatore();
   mappaOrdinata.save("risultati/quadrato1.txt",raw_ascii);
   migliori.save("risultati/miglioriQuadrato.txt",raw_ascii);   
   migliori_semimedia.save("risultati/migliori_semimediaQuadrato.txt",raw_ascii);
   
   return 0;
}

void risolviCommessoViaggiatore(){
	
   cout << "evolvo " << popolazione << " cammini "
	     << " di lunghezza " << Ncities << endl;

   // cammino iniziale
   rowvec cammino = linspace<rowvec>(0,Ncities-1,Ncities);

   // popola la generazione iniziale di cammini
   for (int i=0;i<popolazione;i++) {
      rowvec nuovoCammino = swap(cammino);
      generazione.row(i) = nuovoCammino;
      lunghezze(i) = lunghezzaCammino(nuovoCammino,mappa);
   }
   
   generazione = generazione.rows(sort_index(lunghezze));
	
   for (int k = 0; k<generazioni; k++)  {
   
      // riempio la nuova generazione con i cammini
      // più brevi della vecchia, fatti riprodurre e mutati
      mat nuovaGenerazione = generazione;      
      for (int i=0; i<popolazione/2;i++) if (rnd.Rannyu()<pRiproduzione)
      {
      
         rowvec padre = generazione.row(int(popolazione*
               pow(rnd.Rannyu(),avversitaAmbientale)));
         rowvec madre = generazione.row(int(popolazione*
               pow(rnd.Rannyu(),avversitaAmbientale)));
         mat prole = riproduci( padre, madre );
         nuovaGenerazione.row(i)=muta(prole.row(0));
         nuovaGenerazione.row(i+1)=muta(prole.row(1));
      
      }
      generazione = nuovaGenerazione;
      if( !check(generazione) ) cout <<
         "qualche cammino non soddisfa"<<
         "le condizioni al contorno" << endl;
      
      // ordino la generazione in base alla lunghezza
      for (int i=0; i<popolazione;i++) 
         lunghezze(i) = lunghezzaCammino(generazione.row(i),mappa);
      generazione = generazione.rows(sort_index(lunghezze));
      lunghezze=sort(lunghezze);		
      migliori(k) = int(lunghezze(0));
       migliori_semimedia(k) = int(mean(
         lunghezze.subvec(0, popolazione / 2 - 1)));
         
   }

   // ordino le città in base al cammino minimo 
   // ottenuto come individuo più fit della generazione più evoluta
 
   // mat mappaOrdinata=mappa.rows(sort_index(generazione.row(0)));
   mappaOrdinata=mappa;
   for (int i = 1; i<Ncities; i++)
      mappaOrdinata.row(i)= mappa.row(generazione.row(0)(i)) ;
      
}

//-------------------------------------------------------------
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
      lunghezza += norm(mappa.row(cammino(i))-mappa.row(cammino(i-1)));
   lunghezza += norm(mappa.row(cammino(0))-mappa.row(cammino(Ncities-1)));
   return lunghezza;
   
}

int Pbc(int j) {
	if(j>Ncities-1) return j%Ncities+1;
	return j;
}
 
int randInt(int a, int b){
   return floor(rnd.Rannyu(a,b));
}

void generaMappaCerchio(){
	
	double angle;
	rowvec citta (2);
	for (int i=0;i<Ncities;i++) {
		angle=rnd.Rannyu()*6.28;
		citta(0)=cos(angle);
		citta(1)=sin(angle);
		mappa.row(i)=citta;
	}

	mappa.save("risultati/cerchio0.txt",raw_ascii);
	
}
void generaMappaQuadrato(){

	rowvec citta (2);
	// generate random cities in square 1x1:
	for (int i=0;i<Ncities;i++) {
		citta(0)=rnd.Rannyu();
		citta(1)=rnd.Rannyu();
		mappa.row(i)=citta;
	}
	mappa.save("risultati/quadrato0.txt",raw_ascii);
	
}
//-------------------------------------------------------------
// MUTATIONS: they have to keep the first city in its place
// and the mutated sequence must still be a path

// swap two cities

rowvec muta( rowvec cammino){

         if(rnd.Rannyu()<pScambio) 
            cammino = swap(cammino);  
         if(rnd.Rannyu()<pTrasponi) 
            cammino = trasponi_sottosequenza(cammino);
         if(rnd.Rannyu()<pScambiaSequenze) 
            cammino = scambia_sottosequenze(cammino);
         if(rnd.Rannyu()<pInverti) 
            cammino = inverti_sottosequenza(cammino);
         return cammino;
         
}

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
    return join_vert(camminoPadre, camminoMadre);
}
