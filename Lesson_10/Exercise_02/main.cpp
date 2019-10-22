/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <mpi.h>
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

    MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &Rank);

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

    MPI_Finalize();

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
    //MPI_Barrier(MPI_COMM_WORLD);
    ofstream dist;

    dist.open("output_distance"+to_string(Rank)+".dat",ios::app);

    dist << i << "   " << GetDistance(Path, CityList) << endl;

    if(Rank==0){
        if(i%iprint==0){
            cout<<"step: "<<i<<endl;
            cout<<"temp: "<<temperature<<endl;
        }
    }

    dist.close();
}

void PrintBestPath(vector<int> Path, vector<City> CityList){
    //MPI_Barrier(MPI_COMM_WORLD);

    double dist=GetDistance(Path, CityList);
    double *dist_vec=new double[size];

    MPI_Gather(&dist, 1, MPI_DOUBLE, dist_vec, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    int node;
    if(Rank==0){
        cout<<endl;
        for(int i=0; i<size; i++){
            cout<<"Node "<<i<<" distance: "<<dist_vec[i]<<endl;
        }

        double min_dist=1000;

        for(int i=0; i<size; i++){
            if(dist_vec[i]<min_dist){
                node=i;
                min_dist=dist_vec[i];
            }
        }
        cout<<"Best path found on node "<<node<<endl;
    }

    MPI_Bcast(&node, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if(Rank==node){
        ofstream best_path;
        best_path.open("best_path.dat");

        for(int i=0; i<NUMBER_OF_CITIES; i++){
            int index=Path[i];
            best_path << CityList[index].x << "   " << CityList[index].y <<endl;
        }
        best_path << CityList[Path[0]].x << "   " << CityList[Path[0]].y <<endl;

        best_path.close();
    }

    delete [] dist_vec;

    if(Rank==0)
        rnd.SaveSeed();
}

void Input(){

    int seed[4];
    int * p1;
    int * p2;
    p1=new int[size];
    p2=new int[size];

    if(Rank==0){
        ifstream Primes("Primes");
        if (Primes.is_open()){
            for(int i=0; i<size; i++){
                Primes >> p1[i] >> p2[i] ;
            }
        } else cerr << "PROBLEM: Unable to open Primes" << endl;
        Primes.close();

        ifstream input("seed.in");
        string property;
        if (input.is_open()){
            while ( !input.eof() ){
                input >> property;
                if( property == "RANDOMSEED" ){
                    input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                }
            }
            input.close();
        } else cerr << "PROBLEM: Unable to open seed.in" << endl;

    }

    MPI_Bcast(p1, size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(p2, size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(seed, 4, MPI_INT, 0, MPI_COMM_WORLD);

    rnd.SetRandom(seed, p1[Rank], p2[Rank]);

    if(Rank==0){

        ifstream ReadInput, ReadConfig;

        cout << "The Traveling Salesman Problem  " << endl;
        cout << "Simulated Annealing             " << endl;
        cout << "Parallel MPI code               " << endl << endl;

        //Read input informations
        ReadInput.open("input.dat");

        ReadInput >> T_MAX;
        ReadInput >> T_MIN;
        ReadInput >> T_STEPS;
        ReadInput >> MOVES_PER_TEMPERATURE;

        ReadConfig.open("city_config.dat");
        ReadConfig >> NUMBER_OF_CITIES;

        cout << "Parameters of simulation: " << endl;
        cout << "Number of cities: "<< NUMBER_OF_CITIES << endl;
        cout << "T start: " << T_MAX << endl;
        cout << "T final: " << T_MIN << endl;
        cout << "T steps: " << T_STEPS << endl;
        cout << "Number of moves for each T: " << T_STEPS << endl;
        cout << "Number of nodes requested: " << size << endl << endl;

        ReadInput.close();
        ReadConfig.close();

    }

    MPI_Bcast(&T_MAX, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&T_MIN, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&T_STEPS, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&MOVES_PER_TEMPERATURE, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&NUMBER_OF_CITIES, 1, MPI_INT, 0, MPI_COMM_WORLD);

    coordinate_x=new double[NUMBER_OF_CITIES];
    coordinate_y=new double[NUMBER_OF_CITIES];

    temperature=T_MAX;

    Path.resize(NUMBER_OF_CITIES);
    Path=GeneratePath(NUMBER_OF_CITIES);

    if(Rank==0){
        ifstream ReadConfig;

        ReadConfig.open("city_config.dat");

        ReadConfig >> NUMBER_OF_CITIES; //discard


        for(int i=0; i<NUMBER_OF_CITIES; i++){
            ReadConfig >> coordinate_x[i];
            ReadConfig >> coordinate_y[i];
        }

        ReadConfig.close();
    }

    MPI_Bcast(coordinate_x, NUMBER_OF_CITIES, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(coordinate_y, NUMBER_OF_CITIES, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    //read city positions from file
    Cities.resize(NUMBER_OF_CITIES);
    for(int i=0; i<NUMBER_OF_CITIES; i++){
        Cities[i].x=coordinate_x[i];
        Cities[i].y=coordinate_y[i];
    }

    iprint=T_STEPS/20;

    delete [] coordinate_x;
    delete [] coordinate_y;
    delete [] p1;
    delete [] p2;

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
