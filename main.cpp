
#include "./global.hh"
#include "./ABCDE/parameters.hh"
#include "./ABCDE/WaveVectors.hh"
#include "./ABCDE/omega.hh"
#include "./ABCDE/omega_alphaBN.hh"
#include "./ABCDE/omega_ZnSc.hh"
#include "./ABCDE/omega_alphaBN_bilayer.hh"
#include "./ABCDE/omega_CAC.hh"
#include "./ABCDE/solvediffeq.hh"
#include "./ABCDE/Junc_Density.hh"
#include "./ABCDE/ConcMultiBlock.hh"
#include "./ABCDE/ConcHomo.hh"
#include "./ABCDE/fEhomo.hh"
#include "./ABCDE/Incomp.hh"
#include "./ABCDE/FreeEnergy_Box_Edition.hh"
#include "./ABCDE/size_adjust.hh"
#include "./ABCDE/SaveData.hh"
#include "./ABCDE/FreeEnergy.hh"
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




 //omega_alphaBN(w);
  //omega_CAC(w);
  //omega_alphaBN_bilayer(w);
  //omega_ZnSc(w);

  /*
  double fE_homo;
  std::ofstream outputFile("./RESULTS/p_ave.dat");
  do{
    fE_homo=homogenousfE(chiMatrix,chi);
    std::cout<<Phi_Copo_Dis<<" "<<Phi_Homo_Dis<<" "<<mu_homo<<std::endl;
    outputFile <<Phi_Copo_Dis<<" "<<Phi_Homo_Dis<<" "<<mu_homo<<std::endl;
    mu_homo+=0.1;
  }while(Phi_Homo_Dis<0.999);
  outputFile.close();
  std::cin>>i;
  */
