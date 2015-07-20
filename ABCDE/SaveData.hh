void SaveData(std::vector<double_array> &phi, std::vector<double_array> &w, double_array &dxyz){

  int i, j ,k;

  if(Test==0){
    std::string xyz="./MATLAB/xyz_" + std::to_string(abs(mu_homo)) + "_.dat";
    std::string ABCD="./MATLAB/ABCD_" + std::to_string(abs(mu_homo)) + "_.dat"; 
    //+++++++++++++++++++++++++++++ This output is setup for the matlab plotting +++++++++++++++++++
    std::ofstream outputFile7(xyz);
    for (i=0;i<Nx;i++){
      outputFile7<<i*dxyz(0)<<" "<<i*dxyz(1)<<" "<<i*dxyz(2)<<std::endl;
    }
    outputFile7.close();  
    std::ofstream outputFile8(ABCD);
    outputFile8 <<Phi_Copo_Ord<<" "<<Phi_Homo_Ord<<" "<<0.0<<" "<<0.0<<" "<<0.0<<" "<<0.0<<std::endl;
    for (i=0;i<Nx;i++){
      for(j=0;j<Ny;j++){
	for(k=0;k<Nz;k++){//format A, C, B1+B2+B3, B4
	  outputFile8<<phi[0](i,j,k)<<" "<<phi[1](i,j,k)<<" "<<phi[2](i,j,k)<<" "<<phi[3](i,j,k)<<" "<<phi[4](i,j,k)<<" "<<phi[5](i,j,k)<<std::endl;
	}
      }
    }
    outputFile8.close();
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  }else if(Test == 1){
    std::string xyz="./MATLAB/xyz.dat";
    std::string ABCD="./MATLAB/ABCD.dat"; 
    //+++++++++++++++++++++++++++++ This output is setup for the matlab plotting +++++++++++++++++++
    std::ofstream outputFile7(xyz);
    for (i=0;i<Nx;i++){
      outputFile7<<i*dxyz(0)<<" "<<i*dxyz(1)<<" "<<i*dxyz(2)<<std::endl;
    }
    outputFile7.close();  
    std::ofstream outputFile8(ABCD);
    outputFile8 <<Phi_Copo_Ord<<" "<<Phi_Homo_Ord<<" "<<0.0<<" "<<0.0<<" "<<0.0<<" "<<0.0<<std::endl;
    for (i=0;i<Nx;i++){
      for(j=0;j<Ny;j++){
	for(k=0;k<Nz;k++){//format A, C, B1+B2+B3, B4
	  outputFile8<<phi[0](i,j,k)<<" "<<phi[1](i,j,k)<<" "<<phi[2](i,j,k)<<" "<<phi[3](i,j,k)<<" "<<phi[4](i,j,k)<<" "<<phi[5](i,j,k)<<std::endl;
	}
      }
    }
    outputFile8.close();
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  } else {
    std::cout<<"Somethig is wrong in the SaveData.hh"<<std::endl;
  }
    
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Writting to data files


  std::string omega="./OMEGA/DURING_RUN/omega_" + std::to_string(abs(mu_homo)) + "_.dat";
  std::ofstream outputFile6(omega);
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(k=0;k<Nz;k++){
	outputFile6 <<i<<" "<<j<<" "<<k<< " "<<w[0](i,j,k)<<" "<<w[1](i,j,k)<<" "<<w[2](i,j,k)<<" "<<w[3](i,j,k)<<" "<<w[4](i,j,k)<< " "<<w[5](i,j,k)<<std::endl;
      }
    }
  }
  outputFile6.close();
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Writting Phi1D to data files

  std::string phi_x="./PHI/phix.dat";
  std::string phi_y="./PHI/phiy.dat";
  std::string phi_z="./PHI/phiz.dat";
  
  std::ofstream outputFile16(phi_x);
  for(i=0;i<Nx;i++){
    outputFile16 <<i*dxyz(0)<< " "<<phi[0](i,Ny/2,Nz/2)<<" "<<phi[1](i,Ny/2,Nz/2)<<" "<<phi[2](i,Ny/2,Nz/2)<<" "<<phi[3](i,Ny/2,Nz/2)<<" "<<phi[4](i,Ny/2,Nz/2)<< " "<<phi[5](i,Ny/2,Nz/2)<<std::endl;
  }
  outputFile16.close();

  std::ofstream outputFile26(phi_y);
  for(j=0;j<Ny;j++){
    outputFile26 <<j*dxyz(1)<< " "<<phi[0](Nx/2,j,Nz/2)<<" "<<phi[1](Nx/2,j,Nz/2)<<" "<<phi[2](Nx/2,j,Nz/2)<<" "<<phi[3](Nx/2,j,Nz/2)<<" "<<phi[4](Nx/2,j,Nz/2)<< " "<<phi[5](Nx/2,j,Nz/2)<<std::endl;
  }
  outputFile26.close();
  
  std::ofstream outputFile36(phi_z);
  for(k=0;k<Nz;k++){
    outputFile36 <<k*dxyz(2)<< " "<<phi[0](Nx/2,Ny/2,k)<<" "<<phi[1](Nx/2,Ny/2,k)<<" "<<phi[2](Nx/2,Ny/2,k)<<" "<<phi[3](Nx/2,Ny/2,k)<<" "<<phi[4](Nx/2,Ny/2,k)<< " "<<phi[5](Nx/2,Ny/2,k)<<std::endl;
  }
  outputFile36.close();
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
};



