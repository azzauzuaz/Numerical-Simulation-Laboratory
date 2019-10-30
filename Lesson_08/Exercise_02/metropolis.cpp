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
#include "metropolis.h"

using namespace std;

// NOTE: this code only runs a Monte Carlo simulation for a given set of
//       parameters mu and sigma specified in file input.dat. To run variational
//       optimization of such parameters use Run_variational.sh script instead.

int main (int argc, char *argv[]){

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

    Input();

    for(int iblk=1; iblk <= nblk; ++iblk){ //Simulation

        Reset();

        for(int istep=1; istep <= nstep; ++istep){
            Move();
            if(istep%10 == 0) Measure();
        }

        Averages(iblk);
    }

    rnd.SaveSeed();
    return 0;
}

void Move(){

    xold = x;

    double prob_old=PDF(xold, mu, sigma);

    xnew = x + delta*(rnd.Rannyu() - 0.5) ;

    double prob_new=PDF(xnew, mu, sigma);

    double p = prob_new / prob_old;

    if(p >= rnd.Rannyu()){
        //Update
        x = xnew;
        accepted++;
    }

    attempted++;
}

void Input(void){

    ifstream ReadInput;

    cout << "1D single quantum particle         " << endl;
    cout << "Monte Carlo simulation             " << endl;
    cout << "External potential V(x) = x^4 + 5/2 x^2  " << endl << endl;

    ReadInput.open("input.dat");

    ReadInput >> delta;
    ReadInput >> nblk;
    ReadInput >> nstep;
    ReadInput >> x;
    ReadInput >> mu;
    ReadInput >> sigma;

    hist=new Histogram(-2.5, 2.5, 100, true);

    cout << "The program perform Metropolis moves with uniform translations" << endl;
    cout << "Moves parameter = " << delta << endl;
    cout << "Number of blocks = " << nblk << endl;
    cout << "Number of steps in one block = " << nstep << endl;
    cout << "Starting point of simulation = " << x << endl;
    cout << "mu parameter = " << mu << endl;
    cout << "sigma parameter = " << sigma << endl << endl;

    ReadInput.close();
}

void Reset(){
    sum_h=0;

    attempted = 0;
    accepted = 0;

    blk_norm = 0;
}

void Measure(){
    sum_h+=Hamiltonian(x, mu, sigma);
    hist->add_x(x);
    blk_norm = blk_norm + 1.0;
}

void Averages(int iblk){

    ofstream output_h;

    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << (double)accepted/attempted << endl;

    output_h.open("hamiltonian.dat",ios::app);

    glob_sum_h+=sum_h/blk_norm;
    glob_su2_h+=(sum_h/blk_norm)*(sum_h/blk_norm);

    output_h<< iblk <<"   "<< glob_sum_h/(double)iblk <<"   "<< Error(glob_sum_h,glob_su2_h,iblk) << endl;

    if(iblk==nblk)
        hist->print_hist("sampled_position_histogram.dat");

    output_h.close();
}

double Error(double sum, double sum2, int iblk) {
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}

double PDF(double x, double mu, double sigma){
    double wf=exp(-(x-mu)*(x-mu)/(2.*sigma*sigma))+exp(-(x+mu)*(x+mu)/(2.*sigma*sigma));
    return wf*wf;
}

double Hamiltonian(double x, double mu, double sigma){
    double d2_wf = exp(-(x-mu)*(x-mu)/(2.*sigma*sigma)) * (mu*mu-sigma*sigma+x*x-2.*mu*x)/(sigma*sigma*sigma*sigma) +
                   exp(-(x+mu)*(x+mu)/(2.*sigma*sigma)) * (mu*mu-sigma*sigma+x*x+2.*mu*x)/(sigma*sigma*sigma*sigma);
    double T=x*x*x*x-5./2.*x*x;
    return -0.5 * d2_wf + T;
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
