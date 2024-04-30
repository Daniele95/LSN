
#include "utils.h"

using namespace std;
using namespace arma;


double Norm(vec v) {
   double var = 0.;
   for (int i=0; i<v.size(); i++) var+=v(i)*v(i);
   return sqrt(var);
}

void mediaBlocchi2(vec daMediare, int Nblocchi, 
	int Lstepblocco,string risultati)
{
   vec sum_prog(Nblocchi , arma::fill::zeros);
   vec err_prog(Nblocchi , arma::fill::zeros);
   mediaBlocchi(daMediare, sum_prog,
   	 err_prog, Nblocchi, Lstepblocco); 
   writeMeanError(sum_prog, err_prog,risultati);
}

void writeMeanError(vec sum_prog, vec err_prog, string stringa)
{
   ofstream outfile211(stringa);
   for (int i = 0; i < sum_prog.n_elem; ++i) 
      outfile211 << sum_prog(i) << "\t" << err_prog(i) << endl;
   outfile211.close();
}

void write (float* array, int arrayLength, string nomeFile)
{
   std::ofstream ofile;
   ofile.open(nomeFile);

   for(int i=0; i<arrayLength; i++)
      ofile << array[i] << std::endl;

   ofile.close();   
}

void mediaBlocchi(vec energia, vec& energiaMediaBlocchi, 
	vec& energiaErroreBlocchi, int M_blocchi, int L_dimBlocco) 
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
   	sqrt((mediaQuadra(i) - media(i)*media(i)) / i);
}

float error( float* AV, float* AV2, int n )
{
   if ( n == 0 ) return 0.;
   else return sqrt(  ( AV2[ n ] - pow( AV[ n ], 2) ) / n );
}
                 
