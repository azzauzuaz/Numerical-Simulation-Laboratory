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

double Step(double S, double r, double sigma, double t, double RandomNumber){
    return S*exp((r - 0.5*sigma*sigma)*t + sigma*sqrt(t)*RandomNumber);
};

double CallPayOff(double S, double K){
    double PayOff=S-K;
    if(PayOff>0)
        return PayOff;
    else
        return 0;
}

double PutPayOff(double S, double K){
    double PayOff=K-S;
    if(PayOff>0)
        return PayOff;
    else
        return 0;
}

double error(double AV, double AV2, int n){   // Function for statistical uncertainty estimation
    if(n==0)
        return 0;
    else{
        return sqrt((AV2 - AV*AV)/n);
    }
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

   double S = 100.0;  // asset price at t=0
   double T = 1.0;    // delivery time
   double K = 100.0;  // Strike price
   double r = 0.1;   // risk-free interest rate
   double sigma = 0.25;    // Volatility

   int N_thr=1000;
   int N_blk=100;

   double sum_call, su2_call, sum_put, su2_put;
   //discretized
   double sum_call_disc, su2_call_disc, sum_put_disc, su2_put_disc;

   ofstream output_call("output_call.dat");
   ofstream output_put("output_put.dat");
   ofstream output_call_disc("output_call_discretized.dat");
   ofstream output_put_disc("output_put_discretized.dat");

   for(int i=1; i<N_blk+1; i++){

       sum_call=0;
       su2_call=0;
       sum_put=0;
       su2_put=0;

       sum_call_disc=0;
       su2_call_disc=0;
       sum_put_disc=0;
       su2_put_disc=0;

       for(int j=0; j<N_thr; j++){
           double S_final=Step(S, r, sigma, T, rnd.Gauss(0, 1));
           double S_temp=S;
           for(int k=0; k<100; k++){
               S_temp=Step(S_temp, r, sigma, T/100, rnd.Gauss(0, 1));
           }

           double Call_PO=CallPayOff(S_final, K)*exp(-r*T);
           double Put_PO=PutPayOff(S_final, K)*exp(-r*T);

           double Call_PO_disc=CallPayOff(S_temp, K)*exp(-r*T);
           double Put_PO_disc=PutPayOff(S_temp, K)*exp(-r*T);

           sum_call+=Call_PO;
           su2_call+=Call_PO*Call_PO;
           sum_put+=Put_PO;
           su2_put+=Put_PO*Put_PO;

           sum_call_disc+=Call_PO_disc;
           su2_call_disc+=Call_PO_disc*Call_PO_disc;
           sum_put_disc+=Put_PO_disc;
           su2_put_disc+=Put_PO_disc*Put_PO_disc;

       }

       output_call<<i<<" "<<sum_call/N_thr<<" "<<error(sum_call/N_thr, su2_call/N_thr, i)<<endl;
       output_put<<i<<" "<<sum_put/N_thr<<" "<<error(sum_put/N_thr, su2_put/N_thr, i)<<endl;

       output_call_disc<<i<<" "<<sum_call_disc/N_thr<<" "<<error(sum_call_disc/N_thr, su2_call_disc/N_thr, i)<<endl;
       output_put_disc<<i<<" "<<sum_put_disc/N_thr<<" "<<error(sum_put_disc/N_thr, su2_put_disc/N_thr, i)<<endl;

   }

   output_call.close();
   output_put.close();
   output_call_disc.close();
   output_put_disc.close();

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
