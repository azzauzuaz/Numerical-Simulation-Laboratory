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
#include <vector>
#include "random.h"
#include "variables.h"
#include "functions.h"

using namespace std;

int main (int argc, char *argv[]){

    Input();

    for(int i=0; i<T_STEPS; i++){

        beta=1./temperature;
        Reset();
        for(int j=0; j<MOVES_PER_TEMPERATURE; j++){
            Move();
        }
        PrintDistances(Path, Cities, i);
        temperature=temperature-(T_MAX-T_MIN)/T_STEPS;

    }
    PrintBestPath(Path, Cities);

    return 0;
}

//###############################################################

vector<int> GeneratePath(int NUMBER_OF_CITIES){
    vector<int> Path(NUMBER_OF_CITIES);

    for(int i=0; i<NUMBER_OF_CITIES; i++){
        Path[i] = i;
    }
    for(int i=0; i<NUMBER_OF_CITIES*10; i++){
        Path=Permutation(Path);
    }

    return Path;
}

void Move(){
    vector<int> new_path(Path.size());
    double lold=GetDistance(Path, Cities);

    int rand=rnd.Rannyu(0,5);
    if(rand==0)
        new_path=Permutation(Path);
    if(rand==1)
        new_path=Shift(Path);
    if(rand==2)
        new_path=Partial_Shift(Path);
    if(rand==3)
        new_path=Block_Permutation(Path);
    if(rand==4)
        new_path=Inversion(Path);

    double lnew=GetDistance(new_path, Cities);

    double p = exp(-beta*(lnew-lold));
    
    if(p >= rnd.Rannyu()){
        Path=new_path;
        accepted++;
    }
    attempted++;

}

//mutations
vector<int> Permutation(vector<int> Path){
    vector<int> P_Path=Path;

    int j=rnd.Rannyu(0, Path.size());
    int k=rnd.Rannyu(0, Path.size());

    P_Path[j]=Path[k];
    P_Path[k]=Path[j];

    return P_Path;
}

vector<int> Shift(vector<int> Path){
    int index=rnd.Rannyu(0, Path.size());

    vector<int> Shifted_Path=Path;
    rotate(Shifted_Path.begin(), Shifted_Path.begin()+index, Shifted_Path.end());

    return Shifted_Path;
}

vector<int> Partial_Shift(vector<int> Path){
    int begin=rnd.Rannyu(0, Path.size()/2.);
    int end=rnd.Rannyu(begin, Path.size());

    int index=rnd.Rannyu(begin, end);

    vector<int> Shifted_Path=Path;
    rotate(Shifted_Path.begin()+begin, Shifted_Path.begin()+index ,Shifted_Path.begin()+end);

    return Shifted_Path;
}

vector<int> Block_Permutation(vector<int> Path){
    int half_size=Path.size()/2.;
    int begin=rnd.Rannyu(0, half_size);
    int end=rnd.Rannyu(begin, half_size);

    vector<int> P_Path=Path;

    for(int i=begin; i<end; i++){
        P_Path[i]=Path[i+half_size];
        P_Path[i+half_size]=Path[i];
    }

    return P_Path;
}

vector<int> Inversion(vector<int> Path){
    int begin=rnd.Rannyu(0, Path.size());
    int end=rnd.Rannyu(begin, Path.size());

    vector<int> New_Path=Path;

    for(int i=0; i<end-begin; i++){
        New_Path[i+begin]=Path[end-i-1];
    }

    return New_Path;
}

//fitness
double GetDistance(vector<int> Path, vector<City> CityList){
    double dist=0.;
    int j,k;

    for(int i=0; i<CityList.size()-1; i++){
        j=Path[i];
        k=Path[i+1];
        dist += sqrt( (CityList[j].x-CityList[k].x)*(CityList[j].x-CityList[k].x) +
                      (CityList[j].y-CityList[k].y)*(CityList[j].y-CityList[k].y) );
    }
    j=CityList.size()-1;
    k=0;
    dist += sqrt( (CityList[j].x-CityList[k].x)*(CityList[j].x-CityList[k].x) +
                  (CityList[j].y-CityList[k].y)*(CityList[j].y-CityList[k].y) );

    return dist;
}

double GetDistance2(vector<int> Path, vector<City> CityList){
    double dist=0.;
    int j,k;

    for(int i=0; i<CityList.size()-1; i++){
        j=Path[i];
        k=Path[i+1];
        dist += ( (CityList[j].x-CityList[k].x)*(CityList[j].x-CityList[k].x) +
                  (CityList[j].y-CityList[k].y)*(CityList[j].y-CityList[k].y) );
    }

    j=CityList.size()-1;
    k=0;
    dist += ( (CityList[j].x-CityList[k].x)*(CityList[j].x-CityList[k].x) +
              (CityList[j].y-CityList[k].y)*(CityList[j].y-CityList[k].y) );

    return dist;
}

void PrintDistances(vector<int> Path, vector<City> CityList, int i){
    ofstream dist;

    dist.open("output_distance.dat",ios::app);

    dist << i << "   " << GetDistance(Path, CityList) << endl;

    if(i%iprint==0){
        cout<<"step: "<<i<<endl;
        cout<<"temp: "<<temperature<<endl;
        cout<<"rate: "<<(double)accepted/attempted<<endl<<endl;
    }


    dist.close();
}

void PrintBestPath(vector<int> Path, vector<City> CityList){
    ofstream best_path;
    best_path.open("best_path.dat");

    for(int i=0; i<NUMBER_OF_CITIES; i++){
        int index=Path[i];
        best_path << CityList[index].x << "   " << CityList[index].y <<endl;
    }
    best_path << CityList[Path[0]].x << "   " << CityList[Path[0]].y <<endl;

    best_path.close();

    rnd.SaveSeed();
}

void Input(){

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

    ifstream ReadInput, ReadConfig;

    cout << "The Traveling Salesman Problem  " << endl;
    cout << "Simulated Annealing             " << endl << endl;

    //Read input informations
    ReadInput.open("input.dat");

    ReadInput >> T_MAX;
    ReadInput >> T_MIN;
    ReadInput >> T_STEPS;
    ReadInput >> MOVES_PER_TEMPERATURE;

    temperature=T_MAX;

    ReadConfig.open("city_config.dat");
    ReadConfig >> NUMBER_OF_CITIES;

    Path.resize(NUMBER_OF_CITIES);
    Path=GeneratePath(NUMBER_OF_CITIES);

    //read city positions from file
    Cities.resize(NUMBER_OF_CITIES);
    for(int i=0; i<NUMBER_OF_CITIES; i++){
        ReadConfig >> Cities[i].x;
        ReadConfig >> Cities[i].y;
    }

    cout << "Parameters of simulation: " << endl;
    cout << "Number of cities: "<< NUMBER_OF_CITIES << endl;
    cout << "T start: " << T_MAX << endl;
    cout << "T final: " << T_MIN << endl;
    cout << "T steps: " << T_STEPS << endl;
    cout << "Number of moves for each T: " << T_STEPS << endl<<endl;

    iprint=T_STEPS/20;

    ReadInput.close();
    ReadConfig.close();
}

void Reset(){
    accepted=0;
    attempted=0;
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
