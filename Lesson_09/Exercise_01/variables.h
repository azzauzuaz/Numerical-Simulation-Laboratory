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

struct City{
    double x, y;
};

//random generator
Random rnd;

//variables read from input.dat
int NUMBER_OF_CITIES;
int NUMBER_OF_INDIVIDUALS;
int N_GENERATIONS;
double r;
double p_permutation;
double p_shift;
double p_partial_shift;
double p_block_permutation;
double p_inversion;
double p_crossover;

int iprint;

//vectors for simulation
vector<vector<int> > population;
vector<City> Cities;
vector<vector<int> > new_population;
vector<vector<int> > offspring(2);

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
