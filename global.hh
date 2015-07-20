#include </usr/local/include/fftw3.h>  // This is for My Mac Pro
//#include </opt/sharcnet/fftw/3.3.2/intel/include/fftw3.h> // This is for Sharcnet
//#include </usr/include/fftw3.h> // This is for use on Landua
//#include </usr/local/include/fftw3.h> // This is for elmspace2

#include <stdio.h>     //Include the standard input/output libraries
#include <iostream>  //Cout and Cin etc.
#include <fstream>
#include <stdlib.h>    //Include standard fucntion libraries
#include <math.h>      //Use the math function libraries
#include <time.h>      //Call system time libraries to define the integer seed for random numbers
#include "./include/smemory.hh"  //Use my custom memory handling class
//#include "mpi.h"     //Use this for MPI parallel implimentation later 

#include <vector>

using namespace std;

#define Nx 32
#define Ny 32
#define Nz 32

#define ChainType 6
#define Pi 3.14159

int Iomega;
int box_min;
int Test;
int AlphaB, Bilayer, CAC;


double kappa;
double fA,fC,fB1,fB2,fB3;
double forceD,forceE;
double LAM, HEX, BCC;
double activity;
double mu_homo, mu_copo;
double Phi_Copo_Dis, Phi_Homo_Dis;
double Phi_Copo_Ord, Phi_Homo_Ord;

double Free_Energy, Free_Energy_Homo;
double Lx, Ly, Lz;

double *input_q, *transformed_q, *final_q;
fftw_plan forward_plan, inverse_plan;

typedef array_t<double> double_array;

