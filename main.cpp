#define ARRAY_SIZE(array) (sizeof((array))/sizeof((array[0])))
#include "./global.hh"
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

#include "./HEADERFILES/MatrixInverse.hh"
#include "./HEADERFILES/Inner_Prod_Func.hh"
#include "./HEADERFILES/AndersonMixing.hh"

#include "./HEADERFILES/CalculateFreeEnergy.hh"                        
#include <vector>
#include "./MODS/Mod1.hh"

using namespace std;

// The arguments read in are in order:
// 0-> path (not in use)
// 1-> Type of phase (A=AlphaBN, Z=ZnCs, C=CAC)
// 2-> B2 fraction (Ncopolymer = 200)

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
  DW_0 = create_4d_double_array(History,Nx,Ny,Nz,"DW_0");
  DW_1 = create_4d_double_array(History,Nx,Ny,Nz,"DW_1");
  DW_2 = create_4d_double_array(History,Nx,Ny,Nz,"DW_2");
  DW_3 = create_4d_double_array(History,Nx,Ny,Nz,"DW_3");
  DW_4 = create_4d_double_array(History,Nx,Ny,Nz,"DW_4");
  DW_5 = create_4d_double_array(History,Nx,Ny,Nz,"DW_5");

  W_0 = create_4d_double_array(History,Nx,Ny,Nz,"W_0");
  W_1 = create_4d_double_array(History,Nx,Ny,Nz,"W_1");
  W_2 = create_4d_double_array(History,Nx,Ny,Nz,"W_2");
  W_3 = create_4d_double_array(History,Nx,Ny,Nz,"W_3");
  W_4 = create_4d_double_array(History,Nx,Ny,Nz,"W_4");
  W_5 = create_4d_double_array(History,Nx,Ny,Nz,"W_5");
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
  
  //--------------------------------------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------------------------------------










  
  //--------------------------------------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------------------------------------
  destroy_4d_double_array(DW_0);
  destroy_4d_double_array(DW_1);
  destroy_4d_double_array(DW_2);
  destroy_4d_double_array(DW_3);
  destroy_4d_double_array(DW_4);
  destroy_4d_double_array(DW_5);
  destroy_4d_double_array(W_0);
  destroy_4d_double_array(W_1);
  destroy_4d_double_array(W_2);
  destroy_4d_double_array(W_3);
  destroy_4d_double_array(W_4);
  destroy_4d_double_array(W_5);
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

