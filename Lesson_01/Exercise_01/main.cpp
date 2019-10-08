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

double error(double* AV, double* AV2, int n){   // Function for statistical uncertainty estimation
    if(n==0)
        return 0;
    else{
        //cout<<((AV2[n] - AV[n]*AV[n])/n);
        return sqrt((AV2[n] - AV[n]*AV[n])/n);
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

    // exercise 1.1

    ofstream output1("output01.1.dat");

    int M=100000;              // Total number of throws
    int N=100;                 // Number of blocks
    int L=int(M/N);            // Number of throws in each block, please use for M a multiple of N

    double* r=new double[M];
    for(int i=0; i<M; i++){
        r[i]=rnd.Rannyu();
    }
    double* ave = new double[N];
    double* av2 = new double[N];
    double* sum_prog = new double[N];
    double* su2_prog = new double[N];
    for(int i=0; i<N; i++){
        sum_prog[i]=0;
        su2_prog[i]=0;
    }
    double* err_prog = new double[N];

    for(int i=0; i<N; i++){
        double sum = 0;
        for(int j=0; j<L; j++){
            int k = j+i*L;
            sum += r[k];
        }
        ave[i] = double(sum/L);
        //cout<<ave[i]<<endl;       // r_i
        av2[i] = ave[i]*ave[i]; // (r_i)^2
    }

    for(int i=0; i<N; i++){
        for(int j=0; j<i+1; j++){
            sum_prog[i] += ave[j]; // SUM_{j=0,i} r_j
            su2_prog[i] += av2[j]; // SUM_{j=0,i} (r_j)^2
        }
        sum_prog[i]/=double(i+1); // Cumulative average
        su2_prog[i]/=double(i+1); // Cumulative square average
        err_prog[i] = error(sum_prog,su2_prog,i); // Statistical uncertainty

        output1<<i*L<<"   "<<sum_prog[i]-0.5<<"   "<<err_prog[i]<<endl;

    }

    output1.close();

    // exercise 1.2

    ofstream output2("output01.2.dat");

    for(int i=0; i<M; i++){
        r[i]=rnd.Rannyu();
    }

    for(int i=0; i<N; i++){
        sum_prog[i]=0;
        su2_prog[i]=0;
    }

    for(int i=0; i<N; i++){
        double sum = 0;
        for(int j=0; j<L; j++){
            int k = j+i*L;
            sum += (r[k]-0.5)*(r[k]-0.5);
            //cout<<sum<<endl;
        }
        ave[i] = double(sum/L);
        //cout<<ave[i]<<endl;       // r_i
        av2[i] = ave[i]*ave[i]; // (r_i)^2
    }

    for(int i=0; i<N; i++){
        for(int j=0; j<i+1; j++){
            sum_prog[i] += ave[j]; // SUM_{j=0,i} r_j
            su2_prog[i] += av2[j]; // SUM_{j=0,i} (r_j)^2
        }
        sum_prog[i]/=double(i+1); // Cumulative average
        su2_prog[i]/=double(i+1); // Cumulative square average
        err_prog[i] = error(sum_prog,su2_prog,i); // Statistical uncertainty

        output2<<i*L<<"   "<<sum_prog[i]-1./12.<<"   "<<err_prog[i]<<endl;

    }

    output2.close();

    // exercise 1.3

    ofstream output3("output01.3.dat");

    int m=100;
    int n=1000;
    int jmax=100;

    int* interval=new int[m];

    double chi_avg=0;

    for(int j=0; j<jmax; j++){
        for(int i=0; i<m; i++)
            interval[i]=0;

        for(int i=0; i<n; i++){
            interval[int(rnd.Rannyu()*m)]++;
        }

        double sum=0;
        for(int i=0; i<m; i++){
            sum+=(interval[i]-n/m)*(interval[i]-n/m);
        }

        double chi=sum/(n/m);
        chi_avg+=chi;
        output3<<chi<<endl;
    }

    cout<<chi_avg/m<<endl;

    output3.close();

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
