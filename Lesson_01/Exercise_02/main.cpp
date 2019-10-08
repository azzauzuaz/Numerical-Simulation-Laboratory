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
#include "histogram.h"

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

   Histogram hist1(0, 1, 100, 1);
   Histogram hist2(0, 5, 100, 1);
   Histogram hist3(-10, 10, 100, 1);

   double lambda=1;
   double gamma=1;
   double mu=0;

   int N=100;

   for(int i=0; i<10000; i++){
      double sum1=0;
      double sum2=0;
      double sum3=0;
      for(int j=0; j<N; j++){
          sum1+=rnd.Rannyu();
          sum2+=rnd.Exponential(lambda);
          sum3+=rnd.Cauchy(mu, gamma);
      }

      hist1.add_x(sum1/N);
      hist2.add_x(sum2/N);
      hist3.add_x(sum3/N);
   }

   hist1.print_hist("hist1_out.dat");
   hist2.print_hist("hist2_out.dat");
   hist3.print_hist("hist3_out.dat");

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
