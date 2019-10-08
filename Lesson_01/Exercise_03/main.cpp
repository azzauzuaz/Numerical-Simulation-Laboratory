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
#include <string>
#include <cmath>
#include "random.h"

using namespace std;

int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

   ofstream output("output.dat");

   double L=0.9;
   int Nblk=100;
   int Nthr_per_blk=1000;

   int Nhit;
   double r, x2temp, y2temp, theta;
   double x1, y1, x2, y2;
   double sum=0, su2=0;

   for(int i=1; i<Nblk+1; i++){
       Nhit=0;
       for(int j=0; j<Nthr_per_blk; j++){
           x1=rnd.Rannyu(1, 10);
           y1=rnd.Rannyu(1, 10);

           do{
               x2temp=rnd.Rannyu(-1, 1);
               y2temp=rnd.Rannyu(-1, 1);
               r=x2temp*x2temp+y2temp*y2temp;
           }while(r>1);

           theta=2.*atan(y2temp/(sqrt(x2temp*x2temp+y2temp*y2temp)+x2temp));
           x2=x1+L*cos(theta);
           y2=y1+L*sin(theta);
           if((int)y1!=(int)y2) Nhit++;
       }
       sum+=2.*L*Nthr_per_blk/Nhit;
       su2+=(2.*L*Nthr_per_blk/Nhit)*(2.*L*Nthr_per_blk/Nhit);

       output<<i<<"   "<<sum/i<<"   "<<sqrt( (su2/i - (sum/i)*(sum/i) )/i )<<endl;
   }

   output.close();

   rnd.SaveSeed();
   return 0;
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
