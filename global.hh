#include </usr/local/include/fftw3.h>  // This is for My Mac Pro
//#include </opt/sharcnet/fftw/3.3.2/intel/include/fftw3.h> // This is for Sharcnet
//#include </usr/include/fftw3.h> // This is for use on Landua
//#include </usr/local/include/fftw3.h> // This is for elmspace2

// These were use for the Anderson Mixing section
//#include </System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/Headers/cblas.h>
//#include </System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/Headers/clapack.h>

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "./include/smemory.hh"
#include <vector>
#include <omp.h>

using namespace std;

#define Pi 3.14159

// Setting Nx Ny Nz grid size
#define Nx 32
#define Ny 32
#define Nz 32

// this is used for Andersion Mixing
int Anderson = 0;
int const History = 2;


// Updating parameters
double epsilon_delomega = 0.01;
double epsilon_delphi = 0.01;
double epsilon_delomega_anderson;

// Number of polymer species
#define ChainType 6

// SCFT precision
double  precision=1.0e-4;

int Iomega;
int box_min;
int Test;
int iter;
int AlphaBN, AlphaBNBilayer, Bilayer, CAC, CsCl, ZnSc;
int LAM, HEX, BCC;
int global_index = 0;
double  deltaW;

int NB_middle;
double kappa;
double fA,fC,fB1,fB2,fB3;
double forceD,forceE;

double activity;
double mu_homo, mu_copo;
double Phi_Copo_Dis, Phi_Homo_Dis;
double Phi_Copo_Ord, Phi_Homo_Ord;
double Free_Energy, Free_Energy_Homo;
double Lx, Ly, Lz;
double *input_q, *transformed_q, *final_q;
double *input_q_dag, *transformed_q_dag, *final_q_dag;
int gibberishnamekdshfsjhd = fftw_init_threads();
fftw_plan forward_plan, inverse_plan;
fftw_plan forward_plan_dag, inverse_plan_dag;

typedef array_t<double> double_array;

// Parameters used for Anderson Mixing
double *****DW;
double *****W;
