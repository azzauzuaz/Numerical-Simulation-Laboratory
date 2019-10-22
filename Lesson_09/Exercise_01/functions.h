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

//functions

vector<int> GenerateWorld(int NUMBER_OF_CITIES);
void CheckIndividual(vector<int> Individual);
//selection
int RiggedRoulette(double r, int NUMBER_OF_INDIVIDUALS);
//mutations
vector<int> Permutation(vector<int> Individual);
vector<int> Shift(vector<int> Individual);
vector<int> Partial_Shift(vector<int> Individual);
vector<int> Block_Permutation(vector<int> Individual);
vector<int> Inversion(vector<int> Individual);
vector<vector<int> > Mutation(vector<int> Individual1, vector<int> Individual2);
//crossover
vector<vector<int> > Crossover(vector<int> Individual1, vector<int> Individual2);
//fitness
double GetDistance(vector<int> Individual, vector<City> CityList);
double GetDistance2(vector<int> Individual, vector<City> CityList);
vector<vector<int> > PopulationSort(vector<vector<int> > Individuals, vector<City> CityList);
//utilities
vector<int> sort_indexes(const vector<double> &v);
void Input();
void PrintDistances(vector<vector<int> > Individuals, vector<City> CityList, int generation);
void PrintBestPath(vector<vector<int> > Individuals, vector<City> CityList);

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
