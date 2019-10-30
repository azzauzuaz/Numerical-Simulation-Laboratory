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
#include "metropolis.h"

using namespace std;

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
    cout<<"Equilibration..."<<endl<<endl;
    for(int i=0; i < eqsteps; ++i)
        Move();

    for(int iblk=1; iblk <= nblk; ++iblk){ //Simulation
        Reset();
        for(int istep=1; istep <= nstep; ++istep){
            Move();
            Measure();
        }
        Averages(iblk);
    }

    rnd.SaveSeed();
    return 0;
}

void Move(){

    // 1s move
    xold = x_1s;
    yold = y_1s;
    zold = z_1s;

    double prob_old=PDF_1s(xold,yold,zold);

    //xnew = x_1s + delta_1s*(rnd.Rannyu() - 0.5) ;
    //ynew = y_1s + delta_1s*(rnd.Rannyu() - 0.5) ;
    //znew = z_1s + delta_1s*(rnd.Rannyu() - 0.5) ;

    xnew = x_1s + rnd.Gauss(0, delta_1s);
    ynew = y_1s + rnd.Gauss(0, delta_1s);
    znew = z_1s + rnd.Gauss(0, delta_1s);

    double prob_new=PDF_1s(xnew,ynew,znew);

    double p = prob_new / prob_old;
    if(p >= rnd.Rannyu()){
        //Update
        x_1s = xnew;
        y_1s = ynew;
        z_1s = znew;

        accepted_1s++;
    }
    attempted_1s++;

    //2p move
    xold = x_2p;
    yold = y_2p;
    zold = z_2p;

    prob_old=PDF_2p(xold,yold,zold);

    //xnew = x_2p + delta_2p*(rnd.Rannyu() - 0.5) ;
    //ynew = y_2p + delta_2p*(rnd.Rannyu() - 0.5) ;
    //znew = z_2p + delta_2p*(rnd.Rannyu() - 0.5) ;

    xnew = x_2p + rnd.Gauss(0, delta_2p);
    ynew = y_2p + rnd.Gauss(0, delta_2p);
    znew = z_2p + rnd.Gauss(0, delta_2p);

    prob_new=PDF_2p(xnew,ynew,znew);

    p = prob_new / prob_old;
    if(p >= rnd.Rannyu()){
        //Update
        x_2p = xnew;
        y_2p = ynew;
        z_2p = znew;

        accepted_2p++;
    }
    attempted_2p++;

}

void Input(void){

    double x,y,z;
    ifstream ReadInput;

    cout << "Hydrogen atom                      " << endl;
    cout << "Monte Carlo simulation             " << endl << endl;

    ReadInput.open("input.dat");

    ReadInput >> delta_1s;
    ReadInput >> delta_2p;
    ReadInput >> nblk;
    ReadInput >> nstep;
    ReadInput >> eqsteps;
    ReadInput >> x;
    ReadInput >> y;
    ReadInput >> z;

    x_1s=x;
    y_1s=y;
    z_1s=z;

    x_2p=x;
    y_2p=y;
    z_2p=z;

    cout << "The program perform Metropolis moves with uniform translations" << endl;
    cout << "Moves parameter for 1s= " << delta_1s << endl;
    cout << "Moves parameter for 2p= " << delta_2p << endl;
    cout << "Number of blocks = " << nblk << endl;
    cout << "Number of steps in one block = " << nstep << endl;
    cout << "Starting point of simulation = ( " << x <<", "<< y <<", "<< z <<" )"<< endl << endl;

    ReadInput.close();
}

void Reset(){
    sum_r_1s=0;
    sum_r_2p=0;

    attempted_1s = 0;
    accepted_1s = 0;

    attempted_2p = 0;
    accepted_2p = 0;
}

void Measure(){

    ofstream sampled_pts_1s, sampled_pts_2p;
    ofstream output_position_1s, output_position_2p;

    sampled_pts_1s.open("sampled_1s_points.dat",ios::app);
    sampled_pts_2p.open("sampled_2p_points.dat",ios::app);

    output_position_1s.open("output_1s_position.dat",ios::app);
    output_position_2p.open("output_2p_position.dat",ios::app);

    sampled_pts_1s<<x_1s<<"   "<<y_1s<<"   "<<z_1s<<endl;
    sampled_pts_2p<<x_2p<<"   "<<y_2p<<"   "<<z_2p<<endl;

    double r_1s=sqrt(x_1s*x_1s+y_1s*y_1s+z_1s*z_1s);
    double r_2p=sqrt(x_2p*x_2p+y_2p*y_2p+z_2p*z_2p);

    sum_r_1s+=r_1s;
    sum_r_2p+=r_2p;

    output_position_1s<<r_1s<<endl;
    output_position_2p<<r_2p<<endl;

    sampled_pts_1s.close();
    sampled_pts_2p.close();

    output_position_1s.close();
    output_position_2p.close();
}

void Averages(int iblk){

    ofstream avg_position_1s, avg_position_2p;

    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate for 1s " << (double)accepted_1s/attempted_1s << endl;
    cout << "Acceptance rate for 2p " << (double)accepted_2p/attempted_2p << endl << endl;

    avg_position_1s.open("avg_1s_position.dat",ios::app);
    avg_position_2p.open("avg_2p_position.dat",ios::app);

    glob_sum_r_1s+=sum_r_1s/nstep;
    glob_su2_r_1s+=(sum_r_1s/nstep)*(sum_r_1s/nstep);

    glob_sum_r_2p+=sum_r_2p/nstep;
    glob_su2_r_2p+=(sum_r_2p/nstep)*(sum_r_2p/nstep);

    avg_position_1s<< iblk <<"   "<< glob_sum_r_1s/(double)iblk <<"   "<< Error(glob_sum_r_1s,glob_su2_r_1s,iblk) << endl;
    avg_position_2p<< iblk <<"   "<< glob_sum_r_2p/(double)iblk <<"   "<< Error(glob_sum_r_2p,glob_su2_r_2p,iblk) << endl;

    avg_position_1s.close();
    avg_position_2p.close();
}

double Error(double sum, double sum2, int iblk) {
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}

double PDF_1s(double x, double y, double z){
    double r=sqrt(x*x+y*y+z*z);
    return pow(a0, -3)/M_PI*exp(-2.*r/a0);
}

double PDF_2p(double x, double y, double z){
    double r=sqrt( x*x + y*y + z*z );
    double theta=acos(z/r);
    if(r==0)
        theta=0;
        
    return pow(a0, -5) / (32*M_PI) * r*r * exp(-r/a0)*pow( cos(theta), 2 );
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
