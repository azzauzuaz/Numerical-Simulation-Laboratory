#ifndef __metropolis__
#define __metropolis__

Random rnd;

double x;
double xold;
double xnew;
double mu, sigma;

double sum_h;
double glob_sum_h=0;
double glob_su2_h=0;
double blk_norm;

Histogram* hist;

// simulation
int nstep, nblk;
int accepted, attempted;
double delta;

//functions
void Input(void);
void Move(void);
void Reset(void);
void Measure(void);
void Averages(int);
double Error(double, double, int);
double PDF(double, double, double);
double Hamiltonian(double, double, double);

#endif
