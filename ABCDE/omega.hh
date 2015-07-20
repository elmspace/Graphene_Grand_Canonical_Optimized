void omega(std::vector<double_array> &w){

  int i,j,k;
  int ii,jj,kk;
  double junk;
  double eps =0.001;
 

  if(Iomega==0){
    
    std::ifstream infile;
    if(AlphaB==1){infile.open("./OMEGA/omega_20_20_20_phiC_99.read");}
    if(Bilayer==1){infile.open("./OMEGA/omega_20_20_20_Bilayer.read");}
    if(CAC==1){infile.open("./OMEGA/omega_20_20_20_CAC.read");}

    
    for(i=0;i<Nx;i++){
      for(j=0;j<Ny;j++){
	for(k=0;k<Nz;k++){
	  infile >> ii >> jj >> kk >> w[0](i,j,k) >> w[1](i,j,k) >> w[2](i,j,k) >> w[3](i,j,k) >> w[4](i,j,k) >> w[5](i,j,k);
		
	  if((j>Ny/3)&&(j<2*Ny/3)){
	    
	  }else{
	    w[0](i,j,k)=0.0;
	    w[1](i,j,k)=0.0;
	  }
	}
      }
    }
    infile.close();

  }else if(Iomega==1){
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ Make Structure
    
    for(i=0;i<Nx;i++){
      for(j=0;j<Ny;j++){
	for(k=0;k<Nz;k++){


	  if((i==Nx/2)&&(j==Ny/2)){
	    w[0](i,j,k)=-1.0;//-1.0*cos(k*5.0*Pi/Nz); //A
	    w[1](i,j,k)=-1.0;//1.0*cos(k*5.0*Pi/Nz); //C
	    w[2](i,j,k)=-1.0; //B1
	    w[3](i,j,k)=-1.0; //B2
	    w[4](i,j,k)=-1.0; //B3
	    w[5](i,j,k)=0.0; //B4
	  }else{
	    w[0](i,j,k)=1.0; //A
	    w[1](i,j,k)=1.0; //C
	    w[2](i,j,k)=1.0; //B1
	    w[3](i,j,k)=1.0; //B2
	    w[4](i,j,k)=1.0; //B3
	    w[5](i,j,k)=-10.0; //B4
	  }
	    
	}
      }
    }
 
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  }else if(Iomega==2){
    
    for(i=0;i<Nx;i++){
      for(j=0;j<Ny;j++){
	for(k=0;k<Nz;k++){

	  w[0](i,j,k)=-50.0*(drand48()-0.50); //A1
	  w[1](i,j,k)=-50.0*(drand48()-0.50); //A2
	  w[2](i,j,k)=-50.0*(drand48()-0.50); //B1
	  w[3](i,j,k)=-50.0*(drand48()-0.50); //B2

	}
      }
    }
  


  }
 
};



