#include </usr/local/include/fftw3.h>  // This is for My Mac Pro
//#include </opt/sharcnet/fftw/3.3.2/intel/include/fftw3.h> // This is for Sharcnet
//#include </usr/include/fftw3.h> // This is for use on Landua
//#include </usr/local/include/fftw3.h> // This is for elmspace2

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h> 
#include <math.h> 
#include <time.h>
#include "./include/smemory.hh"
#include <vector>

using namespace std;

#define Pi 3.14159

// Setting Nx Ny Nz grid size
#define Nx 32
#define Ny 32
#define Nz 32

// Updating parameters
double epsilon_delomega = 0.05;
double epsilon_delphi = 0.05;

// Number of polymer species
#define ChainType 6

int Iomega;
int box_min;
int Test;
int AlphaBN, Bilayer, CAC, CsCl, ZnSc;
int LAM, HEX, BCC;

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
fftw_plan forward_plan, inverse_plan;

typedef array_t<double> double_array;

