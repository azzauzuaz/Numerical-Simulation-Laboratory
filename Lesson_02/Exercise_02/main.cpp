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
#include "walker.h"

using namespace std;

double error(double AV, double AV2, int n){   // Function for statistical uncertainty estimation
    if(n==0)
        return 0;
    else
        return sqrt((AV2 - AV*AV)/n);
}

int main (int argc, char *argv[]){

    Discrete_Walker* w1=new Discrete_Walker();
    Continuos_Walker* w2=new Continuos_Walker();

    ofstream output1("output1.dat");
    ofstream output2("output2.dat");

    double AV_1, AV2_1, AV_2, AV2_2;

    for(int step=1; step<=100; step++){

        AV_1=0;
        AV2_1=0;
        AV_2=0;
        AV2_2=0;

        for(int j=0; j<10000; j++){

            w1->reset();
            w2->reset();

            for(int i=0; i<step; i++){
                w1->Walk();
                w2->Walk();
            }

            AV_1+=w1->get_dist();
            AV_2+=w2->get_dist();
            AV2_1+=w1->get_dist()*w1->get_dist();;
            AV2_2+=w2->get_dist()*w2->get_dist();

        }

        AV_1=AV_1/10000;
        AV_2=AV_2/10000;
        AV2_1=AV2_1/10000;
        AV2_2=AV2_2/10000;

        output1<<step<<"    "<<AV_1<<"    "<<error(AV_1, AV2_1 ,step)<<endl;
        output2<<step<<"    "<<AV_2<<"    "<<error(AV_2, AV2_2 ,step)<<endl;
    }

    output1.close();
    output2.close();

    delete w1;
    delete w2;

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
