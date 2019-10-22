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

vector<int> GeneratePath(int NUMBER_OF_CITIES);
//mutations
vector<int> Permutation(vector<int> Path);
vector<int> Shift(vector<int> Path);
vector<int> Partial_Shift(vector<int> Path);
vector<int> Block_Permutation(vector<int> Path);
vector<int> Inversion(vector<int> Path);
//fitness
double GetDistance(vector<int> Path, vector<City> CityList);
double GetDistance2(vector<int> Path, vector<City> CityList);
//utilities
void Input();
void Reset();
void Move();
void PrintDistances(vector<int> Path, vector<City> CityList, int i);
void PrintBestPath(vector<int> Path, vector<City> CityList);

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
