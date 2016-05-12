void omega(std::vector<double_array> &w, double_array &chiMatrix){

  int i,j,k;
  int ii,jj,kk;
  double junk;
  double eps =0.001;
  int fix_Nx=32, fix_Ny=32, fix_Nz=32;
  std::vector<double_array> wT;

  for(int n=0; n<ChainType; n++){
    wT.push_back(double_array(fix_Nx,fix_Ny,fix_Nz));
  }  
  
  
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //++++++++++++++++++++++++++++++++++++   Setting to Zero
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(k=0;k<Nz;k++){
	for(ii=0;ii<ChainType;ii++){
	  w[ii](i,j,k)=0.0;
	}
      }
    }
  }
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++
  


  
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //++++++++++++++++++++++++++++++++++++        Reading In
  if(Iomega==0){

    // ** Read the phi, and not omega, it is more stable
    std::ifstream infile;
    if(AlphaBN==1){infile.open("./OMEGA/READ/omega_AlphaBN_32_32_32.read");}
    if(AlphaBN_single==1){
      if(Nz==16){
	infile.open("./OMEGA/READ/omega_AlphaBN_single_32_32_16.read");
      }else{
	infile.open("./OMEGA/READ/omega_AlphaBN_single_32_32_16.read");
      }
    }
    if(CAC==1){infile.open("./OMEGA/READ/omega_CAC_32_32_32.read");}
    if(CAC_single==1){infile.open("./OMEGA/READ/omega_CAC_single_32_32_32.read");}
    if(ZnSc==1){infile.open("./OMEGA/READ/phi_ZnSc_32_32_32.read");}
    if(ZnSc_single==1){infile.open("./OMEGA/READ/phi_ZnSc_single_32_32_32.read");}
    if(CsCl==1){infile.open("./OMEGA/READ/omega_CsCl_32_32_32.read");}
    if(CsCl_single==1){infile.open("./OMEGA/READ/omega_CsCl_single_32_32_32.read");}
    if(NaCl==1){infile.open("./OMEGA/READ/omega_NaCl_32_32_32.read");}
    if(NaCl_single==1){infile.open("./OMEGA/READ/omega_NaCl_single_32_32_32.read");}
    
    for(i=0;i<fix_Nx;i++){
      for(j=0;j<fix_Ny;j++){
	for(k=0;k<fix_Nz;k+=2){
	  infile >> ii >> jj >> kk >> wT[0](i,j,k) >> wT[1](i,j,k) >> wT[2](i,j,k) >> wT[3](i,j,k) >> wT[4](i,j,k) >> wT[5](i,j,k);
	}
      }
    }
    infile.close();    
  }


  int counter = 1;
  ii=0;
  
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(k=0;k<Nz;k++){

      }
    }
  }


  
 
};



