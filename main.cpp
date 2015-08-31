
#include "./global.hh"
#include "./HEADERFILES/InputArguments.hh"
#include "./HEADERFILES/SetParameters.hh"                             
#include "./HEADERFILES/SetWaveVectors.hh"                             
#include "./HEADERFILES/SetOmega.hh"                                   
#include "./HEADERFILES/SolveModifiedDiffusionEquation.hh"             
#include "./HEADERFILES/CalculateJunctionDensity.hh"                   // <===== ignore for now
#include "./HEADERFILES/ConcMultiBlock.hh"                            
#include "./HEADERFILES/ConcHomo.hh"                                  
#include "./HEADERFILES/CalculateHomogenousFreeEnergy.hh"              
#include "./HEADERFILES/CalculateEta.hh"                              
#include "./HEADERFILES/CalculateFreeEnergy_BoxOptimization.hh"       
#include "./HEADERFILES/OptimizeBoxSize.hh"                            // <===== needs to be checked
#include "./HEADERFILES/SaveData.hh"                                   
#include "./HEADERFILES/CalculateFreeEnergy.hh"                        
#include <vector>
#include "./MODS/Mod1.hh"

using namespace std;

// The arguments read in are in order:
// 0-> path (not in use)
// 1-> Type of phase (A=AlphaBN, Z=ZnCs, C=CAC)
// 2-> B2 fraction (Ncopolymer = 200)

int main(int argc, char* argv[]){


  // --------------------------------------------------------------------------------------------------------------------
  // --------------------------------------------------------------------------------------------------------------------
  // --------------------------------------------------------------------------------------------------------------------
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
  //--------------------------------------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------------------------------------









  //--------------------------------------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------------------------------------
  inputArguments(argc,argv);

  Mod1(w,phi,eta,Ns,ds,k_vector,chi,dxyz,chiMatrix);
  
  //--------------------------------------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------------------------------------










  
  //--------------------------------------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------------------------------------
  fftw_free(input_q);
  fftw_free(transformed_q);
  fftw_free(final_q);
  fftw_destroy_plan(forward_plan);
  fftw_destroy_plan(inverse_plan);
  fftw_cleanup();
  //--------------------------------------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------------------------------------

  return 0;
}

