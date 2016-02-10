#define ARRAY_SIZE(array) (sizeof((array))/sizeof((array[0])))
#include "./global.hh"

// These are the Anderson Mixing (not being used currently)
//#include "./HEADERFILES/MatrixInverse.hh"
//#include "./HEADERFILES/Inner_Prod_Func.hh"
//#include "./HEADERFILES/AndersonMixing.hh"

#include "./HEADERFILES/InputArguments.hh"
#include "./HEADERFILES/SetParameters.hh"
#include "./HEADERFILES/SetWaveVectors.hh"
#include "./HEADERFILES/SetOmega.hh"
#include "./HEADERFILES/SolveModifiedDiffusionEquation_Forward.hh"
#include "./HEADERFILES/SolveModifiedDiffusionEquation_Backward.hh"
#include "./HEADERFILES/CalculateJunctionDensity.hh"
#include "./HEADERFILES/ConcMultiBlock.hh"
#include "./HEADERFILES/ConcHomo.hh"
#include "./HEADERFILES/CalculateHomogenousFreeEnergy.hh"
#include "./HEADERFILES/CalculateEta.hh"
#include "./HEADERFILES/CalculateFreeEnergy_BoxOptimization.hh"
#include "./HEADERFILES/OptimizeBoxSize.hh"
#include "./HEADERFILES/SaveData.hh"
#include "./HEADERFILES/CalculateError.hh"
#include "./HEADERFILES/CalculateAvrg.hh"
#include "./HEADERFILES/SimpleMixing.hh"
#include "./HEADERFILES/CalculateFreeEnergy.hh"
#include <vector>
#include "./MODS/Mod1.hh"
#include "./MODS/Mod0.hh"

using namespace std;

// The arguments read in are in order:
// 0-> path (not in use)
// 1-> Type of phase (look below for the list)
// 2-> B2 fraction (Ncopolymer = 200)


/*
  A -> AlphaBN
  AB -> AlphaBNBilayer
  Z -> ZnSc
  C -> CAC
  L -> Lamellae
  Zc -> ZnSc
  Na -> NaCl
  Cs -> CsCl
*/


int main(int argc, char* argv[]){

  omp_set_num_threads(2);
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


  // Anderson Mixing Parameters *******************************
  create_array(DW,History,ChainType,Nx,Ny,Nz,"DW");
  create_array(W,History,ChainType,Nx,Ny,Nz,"W");
  // ***********************************************************



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

  input_q_dag=(double*)fftw_malloc(sizeof(double)*Nx*Ny*Nz);
  transformed_q_dag=(double*)fftw_malloc(sizeof(double)*Nx*Ny*Nz);
  final_q_dag=(double*)fftw_malloc(sizeof(double)*Nx*Ny*Nz);

  fftw_plan_with_nthreads(2);
  forward_plan=fftw_plan_r2r_3d(Nx,Ny,Nz,input_q,transformed_q,FFTW_REDFT10,FFTW_REDFT10,FFTW_REDFT10,i);
  inverse_plan=fftw_plan_r2r_3d(Nx,Ny,Nz,transformed_q,final_q,FFTW_REDFT01,FFTW_REDFT01,FFTW_REDFT01,i);

  forward_plan_dag=fftw_plan_r2r_3d(Nx,Ny,Nz,input_q_dag,transformed_q_dag,FFTW_REDFT10,FFTW_REDFT10,FFTW_REDFT10,i);
  inverse_plan_dag=fftw_plan_r2r_3d(Nx,Ny,Nz,transformed_q_dag,final_q_dag,FFTW_REDFT01,FFTW_REDFT01,FFTW_REDFT01,i);
  //--------------------------------------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------------------------------------









  //--------------------------------------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------------------------------------
  inputArguments(argc,argv);

  Mod1(w,phi,eta,Ns,ds,k_vector,chi,dxyz,chiMatrix);
  //Mod0(w,phi,eta,Ns,ds,k_vector,chi,dxyz,chiMatrix);

  //--------------------------------------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------------------------------------











  //--------------------------------------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------------------------------------
  destroy_array(W);
  destroy_array(DW);
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
