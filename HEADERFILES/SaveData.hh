void SaveData(std::vector<double_array> &phi, std::vector<double_array> &w, double_array &dxyz){

  int i, j ,k;


  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   Writting Phi for Matlab plotting
  std::string xyz="./MATLAB/xyz.dat";
  std::string ABCD="./MATLAB/ABCD.dat"; 
 
  std::ofstream outputFile1(xyz);
  for (i=0;i<Nx;i++){
    outputFile1<<i*dxyz(0)<<" "<<i*dxyz(1)<<" "<<i*dxyz(2)<<std::endl;
  }
  outputFile1.close();  
  std::ofstream outputFile2(ABCD);
  for (i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(k=0;k<Nz;k++){//format A, C, B1, B2, B3, B4
	outputFile2<<phi[0](i,j,k)<<" "<<phi[1](i,j,k)<<" "<<phi[2](i,j,k)<<" "<<phi[3](i,j,k)<<" "<<phi[4](i,j,k)<<" "<<phi[5](i,j,k)<<std::endl;
      }
    }
  }
  outputFile2.close();
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   Writting omega for read-in later
  std::string omega="./OMEGA/RUN_TIME_DATA/omega.dat";
  std::ofstream outputFile3(omega);
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(k=0;k<Nz;k++){
	outputFile3 <<i<<" "<<j<<" "<<k<< " "<<w[0](i,j,k)<<" "<<w[1](i,j,k)<<" "<<w[2](i,j,k)<<" "<<w[3](i,j,k)<<" "<<w[4](i,j,k)<< " "<<w[5](i,j,k)<<std::endl;
      }
    }
  }
  outputFile3.close();
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   Writting phi for read-in later
  std::string Phi="./OMEGA/RUN_TIME_DATA/phi.dat";
  std::ofstream outputFile4(Phi);
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(k=0;k<Nz;k++){
	outputFile4 <<i<<" "<<j<<" "<<k<< " "<<phi[0](i,j,k)<<" "<<phi[1](i,j,k)<<" "<<phi[2](i,j,k)<<" "<<phi[3](i,j,k)<<" "<<phi[4](i,j,k)<< " "<<phi[5](i,j,k)<<std::endl;
      }
    }
  }
  outputFile4.close();
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   Writting Phi1D to data files
  std::string phi_x="./PHI/phix.dat";
  std::string phi_y="./PHI/phiy.dat";
  std::string phi_z="./PHI/phiz.dat";
  
  std::ofstream outputFile5(phi_x);
  for(i=0;i<Nx;i++){
    outputFile5 <<i*dxyz(0)<< " "<<phi[0](i,Ny/2,Nz/2)<<" "<<phi[1](i,Ny/2,Nz/2)<<" "<<phi[2](i,Ny/2,Nz/2)<<" "<<phi[3](i,Ny/2,Nz/2)<<" "<<phi[4](i,Ny/2,Nz/2)<< " "<<phi[5](i,Ny/2,Nz/2)<<std::endl;
  }
  outputFile5.close();
  
  std::ofstream outputFile6(phi_y);
  for(j=0;j<Ny;j++){
    outputFile6 <<j*dxyz(1)<< " "<<phi[0](Nx/2,j,Nz/2)<<" "<<phi[1](Nx/2,j,Nz/2)<<" "<<phi[2](Nx/2,j,Nz/2)<<" "<<phi[3](Nx/2,j,Nz/2)<<" "<<phi[4](Nx/2,j,Nz/2)<< " "<<phi[5](Nx/2,j,Nz/2)<<std::endl;
  }
  outputFile6.close();
  
  std::ofstream outputFile7(phi_z);
  for(k=0;k<Nz;k++){
    outputFile7 <<k*dxyz(2)<< " "<<phi[0](Nx/2,Ny/2,k)<<" "<<phi[1](Nx/2,Ny/2,k)<<" "<<phi[2](Nx/2,Ny/2,k)<<" "<<phi[3](Nx/2,Ny/2,k)<<" "<<phi[4](Nx/2,Ny/2,k)<< " "<<phi[5](Nx/2,Ny/2,k)<<std::endl;
  }
  outputFile7.close();
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ Writting Phi2D to data files
  std::string phi_xy="./PHI/phixy.dat";
  std::string phi_yz="./PHI/phiyz.dat";
  std::string phi_xz="./PHI/phixz.dat";
  
  std::ofstream outputFile8(phi_xy);
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      outputFile8 <<i*dxyz(0)<<" "<<j*dxyz(1)<< " "<<phi[0](i,j,Nz/2)<<" "<<phi[1](i,j,Nz/2)<<" "<<phi[2](i,j,Nz/2)<<" "<<phi[3](i,j,Nz/2)<<" "<<phi[4](i,j,Nz/2)<< " "<<phi[5](i,j,Nz/2)<<std::endl;
    }
  }
  outputFile8.close();
  std::ofstream outputFile9(phi_yz);
  for(j=0;j<Ny;j++){
    for(k=0;k<Nz;k++){
      outputFile9 <<j*dxyz(1)<<" "<<k*dxyz(2)<< " "<<phi[0](Nx/2,j,k)<<" "<<phi[1](Nx/2,j,k)<<" "<<phi[2](Nx/2,j,k)<<" "<<phi[3](Nx/2,j,k)<<" "<<phi[4](Nx/2,j,k)<< " "<<phi[5](Nx/2,j,k)<<std::endl;
    }
  }
  outputFile9.close();
  std::ofstream outputFile10(phi_xz);
  for(i=0;i<Nx;i++){
    for(k=0;k<Nz;k++){
      outputFile10 <<i*dxyz(0)<<" "<<k*dxyz(2)<< " "<<phi[0](i,Ny/2,k)<<" "<<phi[1](i,Ny/2,k)<<" "<<phi[2](i,Ny/2,k)<<" "<<phi[3](i,Ny/2,k)<<" "<<phi[4](i,Ny/2,k)<< " "<<phi[5](i,Ny/2,k)<<std::endl;
    }
  }
  outputFile10.close();
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

};



