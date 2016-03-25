double ConcMultiBlock(std::vector<double_array> &phi,std::vector<double_array> &w, int *Ns,double ds,double_array &k_vector,double_array &dxyz){

  int         s;
  double      Q;

  double_array	qA(Nx,Ny,Nz,Ns[0]+1);
  double_array	qC(Nx,Ny,Nz,Ns[1]+1);
  double_array	qB1(Nx,Ny,Nz,Ns[2]+1);
  double_array	qB2(Nx,Ny,Nz,Ns[3]+1);
  double_array	qB3(Nx,Ny,Nz,Ns[4]+1);

  double_array	qdagA(Nx,Ny,Nz,Ns[0]+1);
  double_array	qdagC(Nx,Ny,Nz,Ns[1]+1);
  double_array	qdagB1(Nx,Ny,Nz,Ns[2]+1);
  double_array	qdagB2(Nx,Ny,Nz,Ns[3]+1);
  double_array	qdagB3(Nx,Ny,Nz,Ns[4]+1);  

  double_array	qintA(Nx,Ny,Nz);
  double_array	qintC(Nx,Ny,Nz);
  double_array	qintB1(Nx,Ny,Nz);
  double_array	qintB2(Nx,Ny,Nz);
  double_array	qintB3(Nx,Ny,Nz);

  double_array	qintdagA(Nx,Ny,Nz);
  double_array	qintdagC(Nx,Ny,Nz);
  double_array	qintdagB1(Nx,Ny,Nz);
  double_array	qintdagB2(Nx,Ny,Nz);
  double_array	qintdagB3(Nx,Ny,Nz);
 
#pragma omp parallel sections num_threads(2)
  {
    //+++++++++++++++++++++++++++++++++++++++Forward++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#pragma omp section
    {
      qintB1.fill(1.);
      
      solveModDiffEqn_FFT_Forward(qB1,w[2],qintB1,ds,Ns[2],1,k_vector,dxyz);
      for(int i=0;i<Nx;i++){
	for(int j=0;j<Ny;j++){
	  for(int l=0;l<Nz;l++){
	    qintA(i,j,l)=qB1(i,j,l,Ns[2]);
	  }
	}
      }
      
      solveModDiffEqn_FFT_Forward(qA,w[0],qintA,ds,Ns[0],1,k_vector,dxyz);
      for(int i=0;i<Nx;i++){
	for(int j=0;j<Ny;j++){
	  for(int l=0;l<Nz;l++){
	    qintB2(i,j,l)=qA(i,j,l,Ns[0]);
	  }
	}
      }
      
      solveModDiffEqn_FFT_Forward(qB2,w[3],qintB2,ds,Ns[3],1,k_vector,dxyz);
      for(int i=0;i<Nx;i++){
	for(int j=0;j<Ny;j++){
	  for(int l=0;l<Nz;l++){
	    qintC(i,j,l)=qB2(i,j,l,Ns[3]);
	  }
	}
      }
      
      solveModDiffEqn_FFT_Forward(qC,w[1],qintC,ds,Ns[1],1,k_vector,dxyz);
      for(int i=0;i<Nx;i++){
	for(int j=0;j<Ny;j++){
	  for(int l=0;l<Nz;l++){
	    qintB3(i,j,l)=qC(i,j,l,Ns[1]);
	  }
	}
      }
      
      solveModDiffEqn_FFT_Forward(qB3,w[4],qintB3,ds,Ns[4],1,k_vector,dxyz);
    }
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    //+++++++++++++++++++++++++++++++++++++++++++Backward+++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#pragma omp section
    {
      qintdagB3.fill(1.);
      
      solveModDiffEqn_FFT_Backward(qdagB3,w[4],qintdagB3,ds,Ns[4],-1,k_vector,dxyz);
      for(int i=0;i<Nx;i++){
	for(int j=0;j<Ny;j++){
	  for(int l=0;l<Nz;l++){
	    qintdagC(i,j,l)=qdagB3(i,j,l,Ns[4]);
	  }
	}
      }
      
      solveModDiffEqn_FFT_Backward(qdagC,w[1],qintdagC,ds,Ns[1],-1,k_vector,dxyz);
      for(int i=0;i<Nx;i++){
	for(int j=0;j<Ny;j++){
	  for(int l=0;l<Nz;l++){
	    qintdagB2(i,j,l)=qdagC(i,j,l,Ns[1]);
	  }
	}
      }
      
      solveModDiffEqn_FFT_Backward(qdagB2,w[3],qintdagB2,ds,Ns[3],-1,k_vector,dxyz);
      for(int i=0;i<Nx;i++){
	for(int j=0;j<Ny;j++){
	  for(int l=0;l<Nz;l++){
	    qintdagA(i,j,l)=qdagB2(i,j,l,Ns[3]);
	  }
	}
      }
      
      solveModDiffEqn_FFT_Backward(qdagA,w[0],qintdagA,ds,Ns[0],-1,k_vector,dxyz);
      for(int i=0;i<Nx;i++){
	for(int j=0;j<Ny;j++){
	  for(int l=0;l<Nz;l++){
	    qintdagB1(i,j,l)=qdagA(i,j,l,Ns[0]);
	  }
	}
      }
      
      solveModDiffEqn_FFT_Backward(qdagB1,w[2],qintdagB1,ds,Ns[2],-1,k_vector,dxyz);
    }
  }
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  
  //++++++++++++++++++++++++++++++++++++++Single Chain Partition Function+++++++++++++++++++++++++++
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Q=0.0;
#pragma omp parallel for reduction(+:Q)
  for(int i=0;i<Nx;i++){
    for(int j=0;j<Ny;j++){
      for(int l=0;l<Nz;l++){
	Q+=(qdagB1(i,j,l,Ns[2]))*dxyz(0)*dxyz(1)*dxyz(2);
      }
    }
  }
  // Normalizing with respect to the volume of the box
  Q/=((dxyz(0)*Nx)*(dxyz(1)*Ny)*(dxyz(2)*Nz));
  
  
  // Here we do the concentration calculation
  for(int i=0;i<Nx;i++){
    for(int j=0;j<Ny;j++){
      for(int l=0;l<Nz;l++){
	
	
	phi[0](i,j,l)=0.0;  //A
	phi[1](i,j,l)=0.0;  //C
	phi[2](i,j,l)=0.0;  //B1
	phi[3](i,j,l)=0.0;  //B2
	phi[4](i,j,l)=0.0;  //B3
	
	//A1
	for(s=0;s<(Ns[0]+1);s++){
	  if(s==0 || s==Ns[0]){
	    phi[0](i,j,l)+=0.5*qA(i,j,l,s)*qdagA(i,j,l,Ns[0]-s)*ds;
	  }else{
	    phi[0](i,j,l)+=qA(i,j,l,s)*qdagA(i,j,l,Ns[0]-s)*ds;
	  }
	}

	//C
	for(s=0;s<(Ns[1]+1);s++){
	  if(s==0 || s==Ns[1]){
	    phi[1](i,j,l)+=0.5*qC(i,j,l,s)*qdagC(i,j,l,Ns[1]-s)*ds;
	  }else{
	    phi[1](i,j,l)+=qC(i,j,l,s)*qdagC(i,j,l,Ns[1]-s)*ds;
	  }
	}


	//B1
	for(s=0;s<(Ns[2]+1);s++){
	  if(s==0 || s==Ns[2]){
	    phi[2](i,j,l)+=0.5*qB1(i,j,l,s)*qdagB1(i,j,l,Ns[2]-s)*ds;
	  }else{
	    phi[2](i,j,l)+=qB1(i,j,l,s)*qdagB1(i,j,l,Ns[2]-s)*ds;
	  }
	}

	//B2
	for(s=0;s<(Ns[3]+1);s++){
	  if(s==0 || s==Ns[3]){
	    phi[3](i,j,l)+=0.5*qB2(i,j,l,s)*qdagB2(i,j,l,Ns[3]-s)*ds;
	  }else{
	    phi[3](i,j,l)+=qB2(i,j,l,s)*qdagB2(i,j,l,Ns[3]-s)*ds;
	  }
	}

	//B3
	for(s=0;s<(Ns[4]+1);s++){
	  if(s==0 || s==Ns[4]){
	    phi[4](i,j,l)+=0.5*qB3(i,j,l,s)*qdagB3(i,j,l,Ns[4]-s)*ds;
	  }else{
	    phi[4](i,j,l)+=qB3(i,j,l,s)*qdagB3(i,j,l,Ns[4]-s)*ds;
	  }
	}

      }
    }
  }

  if(global_iter%100==0){Junc_Density(qA,qC,qB1,qB2,qB3,qdagA,qdagC,qdagB1,qdagB2,qdagB3,Ns,ds,Q);}
  
  
  

  return Q;


};
