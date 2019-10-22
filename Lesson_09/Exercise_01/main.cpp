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
#include <numeric>      // std::iota
#include <algorithm>    // std::sort
#include "random.h"
#include "variables.h"
#include "functions.h"

using namespace std;

int main (int argc, char *argv[]){

    Input();

    for(int generation=1; generation<=N_GENERATIONS; generation++){
        if(generation%iprint==0)
            cout<<"Generation number: "<<generation<<endl;
        population=PopulationSort(population, Cities);
        for(int j=0; j<NUMBER_OF_INDIVIDUALS; j=j+2){
            int k=RiggedRoulette(r, NUMBER_OF_INDIVIDUALS);
            int l=RiggedRoulette(r, NUMBER_OF_INDIVIDUALS);
            offspring=Crossover(population[k], population[l]);
            offspring=Mutation(offspring[0], offspring[1]);
            CheckIndividual(offspring[0]);
            CheckIndividual(offspring[1]);
            new_population[j]=offspring[0];
            new_population[j+1]=offspring[1];
        }
        population=new_population;
        PrintDistances(population, Cities, generation);

    }
    PrintBestPath(population, Cities);

    return 0;
}

//###############################################################

vector<int> GenerateIndividual(int NUMBER_OF_CITIES){
    vector<int> Individual(NUMBER_OF_CITIES);

    for(int i=0; i<NUMBER_OF_CITIES; i++){
        Individual[i] = i;
    }
    for(int i=0; i<NUMBER_OF_CITIES*10; i++){
        Individual=Permutation(Individual);
    }

    return Individual;
}

void CheckIndividual(vector<int> Individual){
    for(int i=0; i<Individual.size(); i++){
        if(Individual[i]<0 || Individual[i]>Individual.size()){
            cerr<<"Error: individual not fulfilling bonds (wrong value)!"<<endl;
            exit (EXIT_FAILURE);
        }
    }
    for(int i=0; i<Individual.size(); i++){
        for(int j=i+1; j<Individual.size(); j++){
            if(Individual[i]==Individual[j]){
                cerr<<"Error: individual not fulfilling bonds (duplicated value)!"<<endl;
                exit (EXIT_FAILURE);
            }
        }
    }
}

//selection
int RiggedRoulette(double r, int NUMBER_OF_INDIVIDUALS){
    return pow(rnd.Rannyu(), r)*NUMBER_OF_INDIVIDUALS;
}

//mutations
vector<int> Permutation(vector<int> Individual){
    vector<int> P_Individual=Individual;

    int j=rnd.Rannyu(0, Individual.size());
    int k=rnd.Rannyu(0, Individual.size());

    P_Individual[j]=Individual[k];
    P_Individual[k]=Individual[j];

    return P_Individual;
}

vector<int> Shift(vector<int> Individual){
    int index=rnd.Rannyu(0, Individual.size());

    vector<int> Shifted_Individual=Individual;
    rotate(Shifted_Individual.begin(), Shifted_Individual.begin()+index, Shifted_Individual.end());

    return Shifted_Individual;
}

vector<int> Partial_Shift(vector<int> Individual){
    int begin=rnd.Rannyu(0, Individual.size()/2.);
    int end=rnd.Rannyu(begin, Individual.size());

    int index=rnd.Rannyu(begin, end);

    vector<int> Shifted_Individual=Individual;
    rotate(Shifted_Individual.begin()+begin, Shifted_Individual.begin()+index ,Shifted_Individual.begin()+end);

    return Shifted_Individual;
}

vector<int> Block_Permutation(vector<int> Individual){
    int half_size=Individual.size()/2.;
    int begin=rnd.Rannyu(0, half_size);
    int end=rnd.Rannyu(begin, half_size);

    vector<int> P_Individual=Individual;

    for(int i=begin; i<end; i++){
        P_Individual[i]=Individual[i+half_size];
        P_Individual[i+half_size]=Individual[i];
    }

    return P_Individual;
}

vector<int> Inversion(vector<int> Individual){
    int begin=rnd.Rannyu(0, Individual.size());
    int end=rnd.Rannyu(begin, Individual.size());

    vector<int> New_Individual=Individual;

    for(int i=0; i<end-begin; i++){
        New_Individual[i+begin]=Individual[end-i-1];
    }

    return New_Individual;
}

vector<vector<int> > Mutation(vector<int> Individual1, vector<int> Individual2){
    vector<vector<int> > Mutated_Individuals(2);
    Mutated_Individuals[0]=Individual1;
    Mutated_Individuals[1]=Individual2;

    //individual 1
    if(rnd.Rannyu()<p_permutation)
        Mutated_Individuals[0]=Permutation(Mutated_Individuals[0]);
    if(rnd.Rannyu()<p_shift)
        Mutated_Individuals[0]=Shift(Mutated_Individuals[0]);
    if(rnd.Rannyu()<p_partial_shift)
        Mutated_Individuals[0]=Partial_Shift(Mutated_Individuals[0]);
    if(rnd.Rannyu()<p_block_permutation)
        Mutated_Individuals[0]=Block_Permutation(Mutated_Individuals[0]);
    if(rnd.Rannyu()<p_inversion)
        Mutated_Individuals[0]=Inversion(Mutated_Individuals[0]);

    //individual2
    if(rnd.Rannyu()<p_permutation)
        Mutated_Individuals[1]=Permutation(Mutated_Individuals[1]);
    if(rnd.Rannyu()<p_shift)
        Mutated_Individuals[1]=Shift(Mutated_Individuals[1]);
    if(rnd.Rannyu()<p_partial_shift)
        Mutated_Individuals[1]=Partial_Shift(Mutated_Individuals[1]);
    if(rnd.Rannyu()<p_block_permutation)
        Mutated_Individuals[1]=Block_Permutation(Mutated_Individuals[1]);
    if(rnd.Rannyu()<p_inversion)
        Mutated_Individuals[1]=Inversion(Mutated_Individuals[1]);

    return Mutated_Individuals;
}
//fitness
double GetDistance(vector<int> Individual, vector<City> CityList){
    double dist=0.;
    int j,k;

    for(int i=0; i<CityList.size()-1; i++){
        j=Individual[i];
        k=Individual[i+1];
        dist += sqrt( (CityList[j].x-CityList[k].x)*(CityList[j].x-CityList[k].x) +
                      (CityList[j].y-CityList[k].y)*(CityList[j].y-CityList[k].y) );
    }

    j=CityList.size()-1;
    k=0;
    dist += sqrt( (CityList[j].x-CityList[k].x)*(CityList[j].x-CityList[k].x) +
                  (CityList[j].y-CityList[k].y)*(CityList[j].y-CityList[k].y) );

    return dist;
}

double GetDistance2(vector<int> Individual, vector<City> CityList){
    double dist=0.;
    int j,k;

    for(int i=0; i<CityList.size()-1; i++){
        j=Individual[i];
        k=Individual[i+1];
        dist += ( (CityList[j].x-CityList[k].x)*(CityList[j].x-CityList[k].x) +
                  (CityList[j].y-CityList[k].y)*(CityList[j].y-CityList[k].y) );
    }

    j=CityList.size()-1;
    k=0;
    dist += ( (CityList[j].x-CityList[k].x)*(CityList[j].x-CityList[k].x) +
              (CityList[j].y-CityList[k].y)*(CityList[j].y-CityList[k].y) );

    return dist;
}

vector<int> sort_indexes(const vector<double> &v) {

  // initialize original index locations
  vector<int> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

vector<vector<int> > PopulationSort(vector<vector<int> > Individuals, vector<City> CityList){
    vector<double> distances(Individuals.size());

    for(int i=0; i<Individuals.size(); i++){
        distances[i]= GetDistance(Individuals[i], CityList) ;
    }

    vector<int> sorted_indexes=sort_indexes(distances);

    vector<vector<int> > NewPopulation(Individuals.size());

    for(int i=0; i<Individuals.size(); i++){
        NewPopulation[i] = Individuals[sorted_indexes[i]];
    }

    return NewPopulation;
}

//crossover
vector<vector<int> > Crossover(vector<int> Individual1, vector<int> Individual2){
    vector<vector<int> > Offspring(2);
    Offspring[0]=Individual1;
    Offspring[1]=Individual2;

    if(rnd.Rannyu()<p_crossover){
        int index=rnd.Rannyu(0, Individual1.size());

        vector<int> sequence1=Individual2;
        vector<int> sequence2=Individual1;
        for(int i=0; i<sequence1.size(); i++){
            if( find(Individual1.begin()+index, Individual1.end(), sequence1[i]) == Individual1.end() ){//non c'è
                sequence1.erase(sequence1.begin()+i);
                i--;
            }
        }
        for(int i=0; i<sequence2.size(); i++){
            if( find(Individual2.begin()+index, Individual2.end(), sequence2[i]) == Individual2.end() ){//non c'è
                sequence2.erase(sequence2.begin()+i);
                i--;
            }
        }
        for(int i=0; i<sequence1.size(); i++){
            Offspring[0][i+index]=sequence1[i];
            Offspring[1][i+index]=sequence2[i];
        }

    }

    return Offspring;
}

void PrintDistances(vector<vector<int> > Individuals, vector<City> CityList, int generation){
    ofstream best, average;

    best.open("output_best_distance.dat",ios::app);
    average.open("output_average_distance.dat",ios::app);

    vector<vector<int> > sorted_pop=PopulationSort(Individuals, CityList);

    best << generation << "   " << GetDistance(sorted_pop[0], CityList) <<endl;

    double avg_dist=0;
    for(int i=0; i<sorted_pop.size()/2; i++){
        avg_dist+=GetDistance(sorted_pop[i], CityList);
    }

    average << generation << "   " << 2.*avg_dist/sorted_pop.size() <<endl;

    best.close();
    average.close();
}

void PrintBestPath(vector<vector<int> > Individuals, vector<City> CityList){
    ofstream best_path;
    best_path.open("best_path.dat");

    vector<vector<int> > sorted_pop=PopulationSort(Individuals, CityList);

    for(int i=0; i<NUMBER_OF_CITIES; i++){
        int index=sorted_pop[0][i];
        best_path << CityList[index].x << "   " << CityList[index].y <<endl;
    }
    best_path << CityList[sorted_pop[0][0]].x << "   " << CityList[sorted_pop[0][0]].y <<endl;

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
    cout << "Genetic Algorithm               " << endl << endl;

    //Read input informations
    ReadInput.open("input.dat");
    ReadInput >> N_GENERATIONS;
    ReadInput >> NUMBER_OF_INDIVIDUALS;
    ReadInput >> r;
    ReadInput >> p_permutation;
    ReadInput >> p_shift;
    ReadInput >> p_partial_shift;
    ReadInput >> p_block_permutation;
    ReadInput >> p_inversion;
    ReadInput >> p_crossover;

    ReadConfig.open("city_config.dat");
    ReadConfig >> NUMBER_OF_CITIES;

    //populate matrix with individuals on rows and cities on columns
    population.resize(NUMBER_OF_INDIVIDUALS);
    for(int i=0; i<NUMBER_OF_INDIVIDUALS; i++){
        //create and check random individuals
        population[i] = GenerateIndividual(NUMBER_OF_CITIES);
        CheckIndividual(population[i]);
    }

    //read city positions from file
    Cities.resize(NUMBER_OF_CITIES);
    for(int i=0; i<NUMBER_OF_CITIES; i++){
        ReadConfig >> Cities[i].x;
        ReadConfig >> Cities[i].y;
    }

    //new empty population
    new_population.resize(NUMBER_OF_INDIVIDUALS);

    iprint=N_GENERATIONS/10;

    cout << "Parameters of simulation: " << endl;
    cout << "Number of generations: " << N_GENERATIONS << endl;
    cout << "Number of cities: " << NUMBER_OF_CITIES << endl;
    cout << "Population size: " << NUMBER_OF_INDIVIDUALS << endl<<endl;

    ReadInput.close();
    ReadConfig.close();
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
