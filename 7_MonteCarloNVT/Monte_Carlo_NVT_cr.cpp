/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "MD_MC.h"

using namespace std;


//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//parameters, observables
const int m_props=1000; //sovraallocato
int n_props, iv, ik, it, ie, iw, ip;
double vtail, ptail, bin_size, nbins, sd;
double walker[m_props];

// averages
double blk_av[m_props], blk_norm, accepted, attempted;
double glob_av[m_props], glob_av2[m_props];
double stima_pot, stima_pres, stima_kin, stima_etot, stima_temp;
double err_pot, err_pres, err_kin, err_etot, err_temp, err_gdir;

//configuration
const int m_part=108;
double x[m_part],    y[m_part],    z[m_part];
double xold[m_part], yold[m_part], zold[m_part];
double vx[m_part],  vy[m_part],   vz[m_part];

//thermodynamical state
int npart;
double beta,temp,energy,vol,rho,box,rcut;
double Vtail, Wtail;

//simulation
int iNVET, nstep, nblk, restart;
double delta;
std::string phase = " " ;
std::string es;
int ex=0;

//radial distribution function
const int n_dist=100;
double g[n_dist];
double g_ave[n_dist];

//pigreco
const double pi=3.1415927;

//functions
void Input(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
double Boltzmann(double, double, double, int);
double Pbc(double);
double Error(double,double,int);
double Force(int, int);
void Normalization(void);
void Fix_parameters(void);

int main()
{ 

  //Richiesta della fase con cui si vuole fare la simulazione
  int attempt = 0;
  while ( phase != "solid" && phase != "liquid" && phase != "gas" ) {
    if (attempt!=0) { cout << "Phase not valid" << endl; }
    cout << "Insert the phase: solid, liquid or gas: " << endl;
    cin >> phase;
    attempt++;
  }

  //Richiesta dell'esercizio da risolvere
  attempt=0;
  while ( ex != 1 && ex != 2 && ex != 3 && ex != 4 ) {
    if (attempt!=0) { cout << "Exercise not valid" << endl; }
    cout << "Exercise 7.1: digit   --------->   1" << endl; 
    cout << "Exercise 7.2: digit   --------->   2" << endl;
    cout << "Exercise 7.3: digit   --------->   3" << endl;
    cout << "Exercise 7.4: digit   --------->   4" << endl;
    cin >> ex;
    attempt++;
  }

  //Inizialization
  Input(); 

  //Last checkpoint
  bool Test=1;
  while ( Test ) {
    if (ex==4 && iNVET==1) { 
      cout << endl << "This excercise requires Molecular dynamics input, fix it in input."+phase+" file" << endl; 
      return 0;
    }
    else Test=0;

    if (ex!=4 && iNVET==0) { 
      cout << endl << "This excercise requires Monte Carlo input, fix it in input."+phase+" file" << endl; 
      return 0;
    }
    else Test=0;
  }

  //SIMULAZIONE
  Fix_parameters();

  int nconf = 1; 

  //Equilibrazione della temperatura MD NVE
  if(ex==4) {
    for(int iblk=1; iblk <= 20000; iblk++) {
      Reset(iblk);   //Reset block averages
      Move();
    }
  }

  for(int iblk=1; iblk <= nblk; iblk++) {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; istep++) {
      Move();
      Measure();
      Accumulate(); //Update block averages
      if(istep%10 == 0) {
        nconf += 1;
      }
    }
    Averages(iblk);   //Print results for current block
  }

  ConfFinal(); 

  return 0;
}


void Input(void)
{
  ifstream ReadInput, ReadConf, ReadVelocity, Primes, Seed; //Variabili di stream

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "MD(NVE) / MC(NVT) simulation       " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T " << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

//Read seed for random numbers
  int p1, p2;
  Primes.open("Primes");
  Primes >> p1 >> p2 ;
  Primes.close();

//Read input informations
  ReadInput.open("input." + phase);

  ReadInput >> iNVET; //tipologia di simulazione (NVE oppure NVT)
  ReadInput >> restart; //Se restart è 0 (-->Bool=false) allora si riparte dagli stessi numeri di input. Si rifà simulazione identica

  if(restart) Seed.open("seed.out");
  else Seed.open("seed.in");
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed,p1,p2);
  Seed.close();

  ReadInput >> temp; 
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  box = pow(vol,1.0/3.0);
  cout << "Volume of the simulation box = " << vol << endl;
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  cout << "Cutoff of the interatomic potential = " << rcut << endl << endl;
    
  ReadInput >> delta;

  ReadInput >> nblk;

  ReadInput >> nstep;

  cout << "The program perform Metropolis moves with uniform translations" << endl;
  cout << "Moves parameter = " << delta << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();

//Prepare arrays for measurements
  iv = 0; //Potential energy
  it = 1; //Temperature
  ik = 2; //Kinetic energy
  ie = 3; //Total energy
  ip = 4; //Pressure
  n_props = 5; //Number of observables

//Tail corrections
  Vtail = ( (8.0*pi*rho)/(9.0*pow(rcut,9)) - (8.0*pi*rho)/(3.0*pow(rcut,3)) ) * npart;
  Wtail = ( (32.0*pi*rho)/(9.0*pow(rcut,9)) - (16.0*pi*rho)/(3.0*pow(rcut,3)) ) * 3 * npart;

//Read initial configuration
  cout << "Read initial configuration" << endl << endl;
  if(restart) //Se 1 entra nel "se" e parte dalle nuove configurazioni
  {
    ReadConf.open("output/"+phase+"/config_out/config.out");
    ReadVelocity.open("output/"+phase+"/config_out/velocity.out");
    for (int i=0; i<npart; ++i) ReadVelocity >> vx[i] >> vy[i] >> vz[i];
  }
  else 
  {
    ReadConf.open("config.in");
    cout << "Prepare velocities with center of mass velocity equal to zero " << endl;
    double sumv[3] = {0.0, 0.0, 0.0};
    for (int i=0; i<npart; ++i)
    {
      vx[i] = rnd.Gauss(0.,sqrt(temp));
      vy[i] = rnd.Gauss(0.,sqrt(temp));
      vz[i] = rnd.Gauss(0.,sqrt(temp));
      sumv[0] += vx[i];
      sumv[1] += vy[i];
      sumv[2] += vz[i];
    }
    for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart; //Calcola la velocità del centro di massa del sistema
    double sumv2 = 0.0, fs;
    for (int i=0; i<npart; ++i) //Calcola la velocità delle particelle nel sistema del centro di massa
    {
      vx[i] = vx[i] - sumv[0]; 
      vy[i] = vy[i] - sumv[1];
      vz[i] = vz[i] - sumv[2];
      sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
    }
    sumv2 /= (double)npart; 
    fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
    cout << "velocity scale factor: " << fs << endl << endl;
    for (int i=0; i<npart; ++i)
    {
      vx[i] *= fs;
      vy[i] *= fs;
      vz[i] *= fs;
    }
  }

  for (int i=0; i<npart; ++i)
  {
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = Pbc( x[i] * box );
    y[i] = Pbc( y[i] * box );
    z[i] = Pbc( z[i] * box );
  }
  ReadConf.close();

  for (int i=0; i<npart; ++i)
  {
    if(iNVET)
    {
      xold[i] = x[i];
      yold[i] = y[i];
      zold[i] = z[i];
    }
    else
    {
      xold[i] = Pbc(x[i] - vx[i] * delta);
      yold[i] = Pbc(y[i] - vy[i] * delta);
      zold[i] = Pbc(z[i] - vz[i] * delta);
    }
  }
  
//Evaluate properties of the initial configuration
  Measure();

//Print initial values for measured properties
  cout << "Initial potential energy = " << walker[iv]/(double)npart << endl;
  cout << "Initial temperature      = " << walker[it] << endl;
  cout << "Initial kinetic energy   = " << walker[ik]/(double)npart << endl;
  cout << "Initial total energy     = " << walker[ie]/(double)npart << endl;
  cout << "Initial pressure     = " << walker[ip]/(double)npart << endl;

  return;
}


void Move()
{
  int o;
  double p, energy_old, energy_new;
  double xnew, ynew, znew;

  if(iNVET) // Monte Carlo (NVT) move
  {
    for(int i=0; i<npart; ++i)
    {
    //Select randomly a particle (for C++ syntax, 0 <= o <= npart-1)
      o = (int)(rnd.Rannyu()*npart);

    //Old
      energy_old = Boltzmann(x[o],y[o],z[o],o);

    //New
      x[o] = Pbc( x[o] + delta*(rnd.Rannyu() - 0.5) );
      y[o] = Pbc( y[o] + delta*(rnd.Rannyu() - 0.5) );
      z[o] = Pbc( z[o] + delta*(rnd.Rannyu() - 0.5) );

      energy_new = Boltzmann(x[o],y[o],z[o],o);

    //Metropolis test
      p = exp(beta*(energy_old-energy_new));
      if(p >= rnd.Rannyu())  
      {
      //Update
        xold[o] = x[o];
        yold[o] = y[o];
        zold[o] = z[o];
        accepted = accepted + 1.0;
      } else {
        x[o] = xold[o];
        y[o] = yold[o];
        z[o] = zold[o];
      }
      attempted = attempted + 1.0;
    }
  } else // Molecular Dynamics (NVE) move
  {
    double fx[m_part], fy[m_part], fz[m_part];

    for(int i=0; i<npart; ++i){ //Force acting on particle i
      fx[i] = Force(i,0);
      fy[i] = Force(i,1);
      fz[i] = Force(i,2);
    }

    for(int i=0; i<npart; ++i){ //Verlet integration scheme

      xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
      ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
      znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

      vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
      vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
      vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

      xold[i] = x[i];
      yold[i] = y[i];
      zold[i] = z[i];

      x[i] = xnew;
      y[i] = ynew;
      z[i] = znew;

      accepted = accepted + 1.0;
      attempted = attempted + 1.0;
    }
  }
  return;
}

double Boltzmann(double xx, double yy, double zz, int ip)
{
  double ene=0.0;
  double dx, dy, dz, dr;

  for (int i=0; i<npart; ++i)
  {
    if(i != ip)
    {
// distance ip-i in pbc
      dx = Pbc(xx - x[i]);
      dy = Pbc(yy - y[i]);
      dz = Pbc(zz - z[i]);

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);

      if(dr < rcut)
      {
        ene += 1.0/pow(dr,12) - 1.0/pow(dr,6);
      }
    }
  }

  return 4.0*ene;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
 
  return f;
}

void Measure() //Properties measurement
{
  double v = 0.0, kin=0.0, w=0.0;
  double vij;
  double wij; //variabile di appoggio per il calcolo del viriale
  double dx, dy, dz, dr;

  for (int i=0; i<n_dist; i++) { g[i]=0.; }

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i)
  {
    for (int j=i+1; j<npart; ++j)
    {
// distance i-j in pbc
      dx = Pbc(x[i] - x[j]);
      dy = Pbc(y[i] - y[j]);
      dz = Pbc(z[i] - z[j]);

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);

      if(dr < rcut)
      {
        vij = 1.0/pow(dr,12) - 1.0/pow(dr,6);
        v += vij;

        wij = 48 * (1.0/pow(dr,12) - 0.5*1.0/pow(dr,6) );
        w += wij;
      }

      //min_dist = box/2.0;
      int i_bin = floor(2*n_dist*dr/box); //histo position
      if(dr<box/2.){
        g[ i_bin ] += 2;
      }
    }          
  }

  //for (int i=0; i<npart; ++i) kin += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);

  walker[iv] = 4.0 * v + Vtail; // Potential energy
  walker[ip] = w + Wtail ; //calcolo il viriale

  //walker[ik] = kin; // Kinetic energy
  //walker[it] = (2.0 / 3.0) * kin/(double)npart; // Temperature
  //walker[ie] = 4.0 * v + kin;  // Total energy;
  //walker[ip] = rho*walker[it] + 48./(3.*vol)*p/double(npart);


  return;
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
      for(int i=0; i<n_props; ++i)
      {
        glob_av[i] = 0;
        glob_av2[i] = 0;
      }
      for(int i=0; i<n_dist; i++) //Reset radial function
      {
        g[i]=0;
        g_ave[i]=0;
      }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   for(int i=0; i<n_dist; i++) //Reset radial function
   {
     g[i]=0;
     g_ave[i]=0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
   
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i) {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;

   for (int i=0; i<n_dist; i++)
   {
     g_ave[i]=g_ave[i]+g[i];
     //g[i]=0.; //dovrebbe essere già posto a zero in measure()
   }
}


void Averages(int iblk) //Print results for current block
{
    
    ofstream Epot, Pres, G;
    const int wd=20;

    //Nella simulazione con 5e5 blocchi ne stampo uno ogni mille    
    if( es== "/7.2" && iblk%1000 == 0) {
      cout << "----------------------------" << endl << endl;
      cout << "Block number " << iblk << endl;
      cout << "Acceptance rate " << accepted/attempted << endl << endl;
    } 
    if( es!="/7.2" ) {
      cout << "----------------------------" << endl << endl;
      cout << "Block number " << iblk << endl;
      cout << "Acceptance rate " << accepted/attempted << endl << endl;
    }
    
    Epot.open("output/"+phase+es+"/output_epot.dat",ios::app);
    Pres.open("output/"+phase+es+"/output_pres.dat",ios::app);

    if( es!="/7.2" ) {//Calcolo anche la media a blocchi    
      stima_pot = blk_av[iv]/blk_norm/(double)npart ;
      glob_av[iv] += stima_pot;
      glob_av2[iv] += stima_pot*stima_pot;
      err_pot=Error(glob_av[iv],glob_av2[iv],iblk);
    
      stima_pres = rho*temp + (blk_av[ip]/blk_norm)/(3.0*vol); 
      glob_av[ip] += stima_pres;
      glob_av2[ip] += stima_pres*stima_pres;
      err_pres=Error(glob_av[ip],glob_av2[ip],iblk);

      Epot << setw(wd) << iblk <<  setw(wd) << stima_pot << setw(wd) << glob_av[iv]/(double)iblk << setw(wd) << err_pot << endl;
      Pres << setw(wd) << iblk <<  setw(wd) << stima_pres << setw(wd) << glob_av[ip]/(double)iblk << setw(wd) << err_pres << endl;
    }
    else { //Nel punto 7.2 evito di calcolare tutte le medie a blocchi, non servono
      stima_pot = blk_av[iv]/blk_norm/(double)npart ; 
      stima_pres = rho*temp + (blk_av[ip]/blk_norm)/(3.0*vol); 

      Epot << setw(wd) << iblk <<  setw(wd) << stima_pot << endl; 
      Pres << setw(wd) << iblk <<  setw(wd) << stima_pres << endl; 
    }

    Epot.close();
    Pres.close();

    //Radial distribution function
    if (es=="/7.3" || es=="/7.4") {

      G.open("output/"+phase+es+"/output_G.dat",ios::app);
      Normalization();
      
      for(int i=0; i<n_dist; i++) {
        G << g_ave[i] << "  "; //Stampa medie di blocco. Le colonne rappresentano distanze crescenti, le righe rappresentano i diversi blocchi
      }
      G << endl;
      G.close();

    }

}


void ConfFinal(void)
{
  ofstream WriteConf, WriteVelocity, WriteSeed;

  cout << "Print final configuration to file config.out" << endl << endl;
  WriteConf.open("output/"+phase+"/config_out/config.out");
  WriteVelocity.open("output/"+phase+"/config_out/velocity.out");
  for (int i=0; i<npart; ++i)
  {
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
    WriteVelocity << vx[i] << "   " <<  vy[i] << "   " << vz[i] << endl;
  }
  WriteConf.close();
  WriteVelocity.close();

  rnd.SaveSeed();
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r)  //Algorithm for periodic boundary conditions with side L=box
{
    return r - box * rint(r/box);
}

double Error(double sum, double sum2, int iblk)
{
    return sqrt(fabs(sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}

void Normalization() //Normalizza la media a blocchi
{
    double Norm=0;          
    for(int i=0; i<n_dist; i++) {  
      double dr =  box/(2.*n_dist); 
      double r = i*box/(2.*n_dist); 
      Norm = rho*npart * 4.*M_PI/3.*( pow(r+dr,3)-pow(r,3) );
      g_ave[i] = double(g_ave[i]/Norm/(double)blk_norm); 
    }  
}

void Fix_parameters() //Assicura il giusto numero di blocchi anche se ci si dimentica di impostarli
{
    if (ex==1) {
      es="/7.1";
      nstep=1;
      if(phase=="gas") {nblk=2000;} else {nblk=1000;}    
    }

    if (ex==2) {
      es="/7.2";
      nstep=1;
      nblk=5e5;  
    }

    if (ex==3) {
      es="/7.3";
      nstep=2000;
      nblk=50;  
    }

    if (ex==4) {
      es="/7.4";
      nstep=2000;
      nblk=50;  
    }
}
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

