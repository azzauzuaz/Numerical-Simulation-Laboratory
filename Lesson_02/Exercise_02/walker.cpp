#include "walker.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

Walker :: Walker(){
    _x=0;
    _y=0;
    _z=0;

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
};

Walker :: ~Walker(){
      rnd.SaveSeed();
};

void Walker :: reset(){
    _x=0;
    _y=0;
    _z=0;
};

double Walker :: get_dist(){
    return sqrt(_x*_x + _y*_y + _z*_z);

};

void Discrete_Walker::Walk(){
    int pos=rnd.Rannyu(0, 6);
    //cout<<pos<<endl;
    if(pos==0) _x++;
    if(pos==1) _y++;
    if(pos==2) _z++;
    if(pos==3) _x--;
    if(pos==4) _y--;
    if(pos==5) _z--;
};

void Continuos_Walker::Walk(){
    double phi=rnd.Rannyu(0, M_PI);
    double theta=rnd.Rannyu(0, 2.*M_PI);

    _x+=sin(phi)*cos(theta);
    _y+=sin(phi)*sin(theta);
    _z+=cos(phi);

};
