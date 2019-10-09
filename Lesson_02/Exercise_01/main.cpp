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

double fn(double x){
    return M_PI/2.*cos(M_PI*x/2);
}

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

   ofstream output1("output_uniform.dat");
   ofstream output2("output_importance.dat");

   double int_uniform, int_importance;

   int N_blk=100;
   int N_thr=1000;

   double sum_uniform=0;
   double su2_uniform=0;
   double sum_importnace=0;
   double su2_importnace=0;

   for(int i=1; i<N_blk+1; i++){
       int_uniform=0;
       int_importance=0;
       for(int j=0; j<N_thr; j++){
          int_uniform+=fn(rnd.Rannyu());
          double importance_x=1+sqrt(1-rnd.Rannyu());
          int_importance+=fn(importance_x)/(2*(1-importance_x));
       }
       sum_uniform+=int_uniform/N_thr;
       su2_uniform+=(int_uniform/N_thr)*(int_uniform/N_thr);

       sum_importnace+=int_importance/N_thr;
       su2_importnace+=(int_importance/N_thr)*(int_importance/N_thr);

       output1<<i<<"   "<<sum_uniform/i<<"   "<<sqrt( (su2_uniform/i - (sum_uniform/i)*(sum_uniform/i) )/i )<<endl;
       output2<<i<<"   "<<sum_importnace/i<<"   "<<sqrt( (su2_importnace/i - (sum_importnace/i)*(sum_importnace/i) )/i )<<endl;
   }

   output1.close();
   output2.close();

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
