/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

using namespace std;

int size;
int Rank;

struct City{
    double x, y;
};

//random generator
Random rnd;

//variables read from input.dat
int NUMBER_OF_CITIES;
int MOVES_PER_TEMPERATURE;
double T_MIN;
double T_MAX;
int T_STEPS;

double * coordinate_x;
double * coordinate_y;

double beta, temperature;

int accepted, attempted;
int iprint;

//vectors for simulation
vector<int> Path;
vector<City> Cities;

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
