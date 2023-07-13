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
#include <string>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{ 
  Input(); //Inizialization
  while(temp <= ftemp + dtemp * 0.01) //comparing doubles is not really safe...
  {
    ReInput();
    Thermalize(); //Thermalization
    for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
    {
      Reset(iblk);   //Reset block averages
      for(int istep=1; istep <= nstep; ++istep)
      {
        Move(metro);
        Measure();
        Accumulate(); //Update block averages
      }
      Averages(iblk);   //Print results for current block
    }
    ConfFinal(); //Write final configuration
    Results();

    temp += dtemp;
  }

  return 0;
}


void Input(void)
{
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;
  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> ftemp >> dtemp;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> restart; // if=1 restart the simulation from file "config.final"

  ReadInput >> nblk;

  ReadInput >> nstep;

  if(metro==1)
  {
    cout << "The program perform Metropolis moves" << endl;
    path = "metro/";
  }
  else
  {
    cout << "The program perform Gibbs moves" << endl;
    path = "gibbs/";
  }
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();

//Read seed for random numbers
  int p1, p2;
  ifstream Primes("Primes");
  Primes >> p1 >> p2 ;
  Primes.close();

  ifstream Seed;
  if(restart) Seed.open("seed.rest");
  else Seed.open("seed.in");
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed,p1,p2);
  Seed.close();

//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

//initial configuration
  cout << "Read initial configuration" << endl << endl;
  
  if(restart){
    ifstream Conf("config.rest");
    for (int i=0; i<nspin; ++i) Conf >> s[i];
    Conf.close();
  }
  else
  {
    for (int i=0; i<nspin; ++i)
    {
      if(rnd.Rannyu() >= 0.5) s[i] = 1;
      else s[i] = -1;
    }
  }
  
//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl << endl;
}


void ReInput(void)
{
  cout << "============================" << endl << endl;
  cout << "Current temperature = " << temp << endl << endl;

  beta = 1.0/temp;

  for (int i=0; i<nspin; ++i)
  {
    if(rnd.Rannyu() >= 0.5) s[i] = 1;
    else s[i] = -1;
  }
}


void Thermalize(){
  const int Nterm=2000;
  cout << "Thermalization" << endl;
  for(int i=0; i<Nterm; ++i) Move(metro);
  cout << "Done" << endl << endl;
}


void Move(int metro)
{
  int o;
  double p, energy_old, energy_new, sm;
  double energy_up, energy_down;

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);

    if(metro==1) //Metropolis
    {
      sm = s[o];

    //Old
      energy_old = Boltzmann(sm,o);

    //New
      sm *= -1;
      energy_new = Boltzmann(sm,o);

    //Metropolis test
      p = exp(beta*(energy_old-energy_new)); //p>1 if energy_new < energy_old
      if(rnd.Rannyu() <= p)  
      {
      //Update
        s[o] = sm;
        accepted++;
      }
      attempted++;
    }
    else //Gibbs sampling
    {
      energy_up = Boltzmann(1,o);
      energy_down = Boltzmann(-1,o);
      p = 1. / (1 + exp(beta*(energy_up-energy_down)));

      sm = -1;
      if (rnd.Rannyu() <= p) sm = 1;
      s[o] = sm;

      accepted++;
      attempted++;
    }
  }
}


double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}


void Measure()
{
  //int bin;
  double u = 0.0, m = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
     m += s[i];
  }
  walker[iu] = u;
  walker[ic] = u * u; //at the end of each block we will compute beta^2(<u^2>-<u>^2)
  walker[im] = m;
  walker[ix] = m * m; //multiply by beta at the end
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
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    
   ofstream Ene, Heat, Mag, Chi;
   const int wd=13;
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    Ene.open(path+"output.ene.0",ios::app);
    stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    err[iu]=Error(glob_av[iu],glob_av2[iu],iblk);
    Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err[iu] << endl;
    Ene.close();

    Heat.open(path+"output.heat.0",ios::app);
    stima_c = beta*beta*(blk_av[ic] - blk_av[iu]*blk_av[iu]/blk_norm)/blk_norm/(double)nspin; //Heat Capacity
    glob_av[ic]  += stima_c;
    glob_av2[ic] += stima_c*stima_c;
    err[ic]=Error(glob_av[ic],glob_av2[ic],iblk);
    Heat << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err[ic] << endl;
    Heat.close();

    Mag.open(path+"output.mag.0",ios::app);
    stima_m = blk_av[im]/blk_norm/(double)nspin; //Magnetization
    glob_av[im]  += stima_m;
    glob_av2[im] += stima_m*stima_m;
    err[im]=Error(glob_av[im],glob_av2[im],iblk);
    Mag << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err[im] << endl;
    Mag.close();

    Chi.open(path+"output.chi.0",ios::app);
    stima_x = beta*blk_av[ix]/blk_norm/(double)nspin; //Susceptibility
    glob_av[ix]  += stima_x;
    glob_av2[ix] += stima_x*stima_x;
    err[ix]=Error(glob_av[ix],glob_av2[ix],iblk);
    Chi << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err[ix] << endl;
    Chi.close();

    cout << "----------------------------" << endl << endl;
}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}


void Results(void)
{
  ofstream Ene, Heat, Mag, Chi;
  const int wd=13;

  cout << "Print final values" << endl << endl;

  if (h == 0)
  {
    Ene.open(path+"output.ene.final",ios::app);
    Ene << setw(wd) << temp << setw(wd) << glob_av[iu]/(double)nblk << setw(wd) << err[iu] << endl;
    Ene.close();

    Heat.open(path+"output.heat.final",ios::app);
    Heat << setw(wd) << temp << setw(wd) << glob_av[ic]/(double)nblk << setw(wd) << err[ic] << endl;
    Heat.close();

    Chi.open(path+"output.chi.final",ios::app);
    Chi << setw(wd) << temp << setw(wd) << glob_av[ix]/(double)nblk << setw(wd) << err[ix] << endl;
    Chi.close();
  }
  else
  {
  Mag.open(path+"output.mag.final",ios::app);
  Mag << setw(wd) << temp << setw(wd) << glob_av[im]/(double)nblk << setw(wd) << err[im] << endl;
  Mag.close();
  }

}


int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}


double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
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
