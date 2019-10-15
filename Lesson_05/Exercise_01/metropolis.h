#ifndef __metropolis__
#define __metropolis__

Random rnd;

//const double a0 = 5.29E-11;
const double a0 = 1; //Use Bohr radius units, a0 for distances

double x_1s, x_2p;
double y_1s, y_2p;
double z_1s, z_2p;

double xold;
double yold;
double zold;

double xnew;
double ynew;
double znew;

double sum_r_1s, sum_r_2p;
double glob_sum_r_1s=0;
double glob_su2_r_1s=0;
double glob_sum_r_2p=0;
double glob_su2_r_2p=0;

// simulation
int nstep, nblk, eqsteps;
int accepted_1s, attempted_1s, accepted_2p, attempted_2p;
double delta_1s, delta_2p;

//functions
void Input(void);
void Move(void);
void Reset(void);
void Measure(void);
void Averages(int);
double Error(double, double, int);
double PDF_1s(double, double, double);
double PDF_2p(double, double, double);

#endif
