#include "./OMEGA/SetOmega_alphaBN_bilayer.hh"
#include "./OMEGA/SetOmega_alphaBN.hh"
#include "./OMEGA/SetOmega_CAC.hh"
#include "./OMEGA/SetOmega_ZnSc.hh"


void omega(std::vector<double_array> &w, double_array &chiMatrix){

  int i,j,k;
  int ii,jj,kk;
  double junk;
  double eps =0.001;
  std::vector<double_array> phi_w;
  for(int n=0; n<ChainType; n++) 
    {
      phi_w.push_back(double_array(Nx,Ny,Nz));
    }  

  if(Iomega==0){

    // ** Read the phi, and not omega, it is more stable
    std::ifstream infile;
    if(AlphaB==1){infile.open("./OMEGA/omega_20_20_20_phiC_99.read");}
    if(Bilayer==1){infile.open("./OMEGA/omega_20_20_20_Bilayer.read");}
    if(CAC==1){infile.open("./OMEGA/omega_20_20_20_CAC.read");}

    for(i=0;i<Nx;i++){
      for(j=0;j<Ny;j++){
	for(k=0;k<Nz;k++){
	  infile >> ii >> jj >> kk >> phi_w[0](i,j,k) >> phi_w[1](i,j,k) >> phi_w[2](i,j,k) >> phi_w[3](i,j,k) >> phi_w[4](i,j,k) >> phi_w[5](i,j,k);	  
	}
      }
    }
    infile.close();

    // Setting the omega field from the input phi
    for(i=0;i<Nx;i++){
      for(j=0;j<Ny;j++){
	for(k=0;k<Nz;k++){
	  for(ii=0;ii<ChainType;ii++){
	    w[ii](i,j,k)=0.0;
	    for(jj=0;jj<ChainType;jj++){
	      w[ii](i,j,k)+=chiMatrix(ii,jj)*phi_w[jj](i,j,k);
	    }
	  }
	}
      }
    }
    

  }else if(Iomega==1){
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ Make Structure
    // Calling a given function for setting up the phi and omega fields
    /*
    if(alphaBN==1){}
    if(alphaBN_Bilayer==1){}
    if(CAC==1){}
    if(LAM==1){}
    if(HEX==1){}
    if(BCC==1){}   
    */
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  }







  
 
};



