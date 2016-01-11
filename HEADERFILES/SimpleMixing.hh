/*
  Here we have the Simple Mixing algorithm
*/
using namespace std;
void SimpleMixing(std::vector<double_array> &w, std::vector<double_array> &newW, std::vector<double_array> &delW,double_array &delphi, double_array &dxyz){

  int h, i, j, k, s;
  int chain;
  
  // If Anderson Mixing is turned off, then simple mixing is used
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(k=0;k<Nz;k++){
	for(chain=0;chain<ChainType;chain++){
	  w[chain](i,j,k)+=(epsilon_delomega*delW[chain](i,j,k)-epsilon_delphi*delphi(i,j,k));
	}
      }
    }
  }
  
}
