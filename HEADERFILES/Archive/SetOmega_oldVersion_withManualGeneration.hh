void omega(std::vector<double_array> &w, double_array &chiMatrix){

  int i,j,k;
  int ii,jj,kk;
  double junk;
  double eps =0.001;
  std::vector<double_array> phi_w;

  for(int n=0; n<ChainType; n++){
    phi_w.push_back(double_array(Nx,Ny,Nz));
  }  
  
  
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //++++++++++++++++++++++++++++++++++++   Setting to Zero
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(k=0;k<Nz;k++){
	for(ii=0;ii<ChainType;ii++){
	  phi_w[ii](i,j,k)=0.0;
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
	infile.open("./OMEGA/READ/omega_AlphaBN_single_32_32_32.read");
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
    
    for(i=0;i<Nx;i+=1){
      for(j=0;j<Ny;j+=1){
	for(k=0;k<Nz;k+=1){
	  //infile >> ii >> jj >> kk >> phi_w[0](i,j,k) >> phi_w[1](i,j,k) >> phi_w[2](i,j,k) >> phi_w[3](i,j,k) >> phi_w[4](i,j,k) >> phi_w[5](i,j,k);
	  infile >> ii >> jj >> kk >> w[0](i,j,k) >> w[1](i,j,k) >> w[2](i,j,k) >> w[3](i,j,k) >> w[4](i,j,k) >> w[5](i,j,k);
	  if(AlphaBN==1){w[5](i,j,k)=0.0;}
	}
      }
    }
    infile.close();    

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++         Setting W
  }else if(Iomega==1){
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++               LAM
    if(LAM==1){

      for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
	  for(k=0;k<Nz;k++){
	    phi_w[0](i,j,k)=cos(2.0*Pi*i/Nx)+1.0;
	  }
	}
      }
        
    }
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++               HEX
    if(HEX==1){

      for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
	  for(k=0;k<Nz;k++){
	    phi_w[0](i,j,k)=(cos(2.0*Pi*i/Nx)*cos(2.0*Pi*j/Ny))+1.0; 
	  }
	}
      }
        
    }
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++               BCC
    if(BCC==1){

      phi_w[0](0,0,0)=1.0;
      phi_w[0](Nx-1,0,0)=1.0;
      phi_w[0](0,Ny-1,0)=1.0;
      phi_w[0](0,0,Nz-1)=1.0;
      phi_w[0](Nx-1,Ny-1,0)=1.0;
      phi_w[0](0,Ny-1,Nz-1)=1.0;
      phi_w[0](Nx-1,0,Nz-1)=1.0;
      phi_w[0](Nx-1,Ny-1,Nz-1)=1.0;
      phi_w[0]((Nx-1)/2,(Ny-1)/2,(Nz-1)/2)=1.0;
      
        
    }
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++              CsCl
    if(BCC==1){

      phi_w[0](0,0,0)=1.0;
      phi_w[0](Nx-1,0,0)=1.0;
      phi_w[0](0,Ny-1,0)=1.0;
      phi_w[0](0,0,Nz-1)=1.0;
      phi_w[0](Nx-1,Ny-1,0)=1.0;
      phi_w[0](0,Ny-1,Nz-1)=1.0;
      phi_w[0](Nx-1,0,Nz-1)=1.0;
      phi_w[0](Nx-1,Ny-1,Nz-1)=1.0;
      phi_w[1]((Nx-1)/2,(Ny-1)/2,(Nz-1)/2)=1.0;
      
    }
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++              ZnSc
    if(ZnSc==1){

      /*
      // Corners
      phi_w[0](0,0,0)=100.0;
      phi_w[0](Nx-1,0,0)=100.0;
      phi_w[0](0,Ny-1,0)=100.0;
      phi_w[0](0,0,Nz-1)=100.0;
      phi_w[0](Nx-1,Ny-1,0)=100.0;
      phi_w[0](0,Ny-1,Nz-1)=100.0;
      phi_w[0](Nx-1,0,Nz-1)=100.0;
      phi_w[0](Nx-1,Ny-1,Nz-1)=100.0;
      // Faces
      phi_w[0](0,Ny/2,Nz/2)=100.0;
      phi_w[0](Nx/2,0,Nz/2)=100.0;
      phi_w[0](Nx/2,Ny/2,0)=100.0;
      phi_w[0](Nx-1,Ny/2,Nz/2)=100.0;
      phi_w[0](Nx/2,Ny-1,Nz/2)=100.0;
      phi_w[0](Nx/2,Ny/2,Nz-1)=100.0;

      phi_w[1](Nx/4,Ny/4,3*Nz/4)=1000.0;
      phi_w[1](Nx/4,3*Ny/4,Nz/4)=1000.0;
      phi_w[1](3*Nx/4,Ny/4,Nz/4)=1000.0;
      phi_w[1](3*Nx/4,3*Ny/4,3*Nz/4)=1000.0;
      */

      junk=-20.0;
      // Corners
      w[0](0,0,0)=junk;
      w[0](Nx-1,0,0)=junk;
      w[0](0,Ny-1,0)=junk;
      w[0](0,0,Nz-1)=junk;
      w[0](Nx-1,Ny-1,0)=junk;
      w[0](0,Ny-1,Nz-1)=junk;
      w[0](Nx-1,0,Nz-1)=junk;
      w[0](Nx-1,Ny-1,Nz-1)=junk;
      // Faces
      w[0](0,Ny/2,Nz/2)=junk;
      w[0](Nx/2,0,Nz/2)=junk;
      w[0](Nx/2,Ny/2,0)=junk;
      w[0](Nx-1,Ny/2,Nz/2)=junk;
      w[0](Nx/2,Ny-1,Nz/2)=junk;
      w[0](Nx/2,Ny/2,Nz-1)=junk;

      w[1](Nx/4,Ny/4,3*Nz/4)=junk;
      w[1](Nx/4,3*Ny/4,Nz/4)=junk;
      w[1](3*Nx/4,Ny/4,Nz/4)=junk;
      w[1](3*Nx/4,3*Ny/4,3*Nz/4)=junk;
      
        
    }
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++               CAC
    if(CAC==1){
      
      for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
	  for(k=0;k<Nz;k++){
	    phi_w[0](i,j,k)=(cos(2.0*Pi*i/Nx)*cos(2.0*Pi*j/Ny))+1.0;
	    phi_w[1](i,j,k)=(cos(2.0*Pi*i/Nx)*cos(2.0*Pi*j/Ny))+1.0; 
	  }
	}
      }
        
    }
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++


    //++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++           AlphaBN
    if(AlphaBN==1){

      for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
	  for(k=0;k<4;k++){
	    phi_w[0](i,j,k)=(cos(2.0*Pi*i/Nx)*cos(2.0*Pi*j/Ny))+1.0; 
	  }
	}
      }

      for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
	  for(k=(Nz-4);k<Nz;k++){
	    phi_w[0](i,j,k)=(cos(2.0*Pi*i/Nx)*cos(2.0*Pi*j/Ny))+1.0; 
	  }
	}
      }
 
      for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
	  for(k=(Nz/2-2);k<(Nz/2+3);k++){
	    phi_w[1](i,j,k)=(cos(2.0*Pi*i/Nx)*cos(2.0*Pi*j/Ny))+1.0;
	    phi_w[5](i,j,k)=-1.0;
	  }
	}
      }

      for(k=(Nz/2-2);k<(Nz/2+3);k++){
	phi_w[1](Nx/4,0,0)=5.0;
	phi_w[1](Nx/4,Ny-1,0)=5.0;
	phi_w[1](3*Nx/4,Ny/2,0)=5.0;
	
	phi_w[1](Nx/4,0,Nz-1)=5.0;
	phi_w[1](Nx/4,Ny-1,Nz-1)=5.0;
	phi_w[1](3*Nx/4,Ny/2,Nz-1)=5.0;
	
	phi_w[0](Nx/4,0,k)=5.0;
	phi_w[0](Nx/4,Ny-1,k)=5.0;
	phi_w[0](3*Nx/4,Ny/2,k)=5.0;
      }

      
    }
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++


  }


  /*
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */

  
 
};



