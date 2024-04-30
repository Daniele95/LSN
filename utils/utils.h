#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <armadillo>

using namespace std;
using namespace arma;

double Norm(vec v) ;
void mediaBlocchi(vec energia, vec& energiaMediaBlocchi, 
	vec& energiaErroreBlocchi, int M_blocchi, int L_dimBlocco); 
 
double error(vec media, vec mediaQuadra, int i) ;

float error( float* AV, float* AV2, int n );
void writeMeanError(vec sum_prog, 
	vec err_prog, string stringa);
void mediaBlocchi2(vec daMediare, int Nblocchi, 
	int Lstepblocco,string risultati);
	
void write (float* array, int arrayLength, string nomeFile);
#endif
