
#include "./global.hh"
#include "./HEADERFILES/SetParameters.hh"                              // <===== needs fixing
#include "./HEADERFILES/SetWaveVectors.hh"                             // <===== good
#include "./HEADERFILES/SetOmega.hh"                                   // <===== needs fixing
#include "./HEADERFILES/SolveModifiedDiffusionEquation.hh"             // <===== good
#include "./HEADERFILES/CalculateJunctionDensity.hh"                   // <===== ignore for now
#include "./HEADERFILES/ConcMultiBlock.hh"                             // <===== good
#include "./HEADERFILES/ConcHomo.hh"                                   // <===== good
#include "./HEADERFILES/CalculateHomogenousFreeEnergy.hh"              // <===== good
#include "./HEADERFILES/CalculateEta.hh"                               // <===== good
#include "./HEADERFILES/CalculateFreeEnergy_BoxOptimization.hh"        // <===== good
#include "./HEADERFILES/OptimizeBoxSize.hh"                            // <===== needs to be checked
#include "./HEADERFILES/SaveData.hh"                                   // <===== good
#include "./HEADERFILES/CalculateFreeEnergy.hh"                        // <===== good
#include <vector>
#include "./MODS/Mod1test.hh"

using namespace std;


int main(){
  std::vector<double_array> w;
  std::vector<double_array> phi;
  double_array eta(Nx,Ny,Nz);
  double_array chi(15);	// There are only 15 interactions
  double_array f(ChainType);
  int * Ns;
  double_array k_vector(Nx,Ny,Nz);
  double_array dxyz(3);
  double_array chiMatrix(ChainType,ChainType);
  w.reserve(ChainType);
  phi.reserve(ChainType);

  for(int n=0; n<ChainType; n++) 
  {
	w.push_back(double_array(Nx,Ny,Nz));
	phi.push_back(double_array(Nx,Ny,Nz));
  }  
  
  Ns=new int[ChainType];

  int i;
  double  ds;

  long iseed;
  time_t t;
  iseed=time(&t);
  
  srand48(iseed);

  input_q=(double*)fftw_malloc(sizeof(double)*Nx*Ny*Nz);
  transformed_q=(double*)fftw_malloc(sizeof(double)*Nx*Ny*Nz);
  final_q=(double*)fftw_malloc(sizeof(double)*Nx*Ny*Nz);


  forward_plan=fftw_plan_r2r_3d(Nx,Ny,Nz,input_q,transformed_q,FFTW_REDFT10,FFTW_REDFT10,FFTW_REDFT10,i);
  inverse_plan=fftw_plan_r2r_3d(Nx,Ny,Nz,transformed_q,final_q,FFTW_REDFT01,FFTW_REDFT01,FFTW_REDFT01,i);
  
 
  //________________________________________________________________________________________________________________________________

  Mod1(w,phi,eta,Ns,ds,k_vector,chi,dxyz,chiMatrix);
  
  //________________________________________________________________________________________________________________________________
  

  //Destroy memory allocations------------
  fftw_free(input_q);
  fftw_free(transformed_q);
  fftw_free(final_q);
  
  fftw_destroy_plan(forward_plan);
  fftw_destroy_plan(inverse_plan);
  fftw_cleanup();

  //-------------------------------------


  return 0;
}

