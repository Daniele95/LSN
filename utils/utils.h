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


#endif
