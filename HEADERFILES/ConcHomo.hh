double ConcHomo(std::vector<double_array> &phi,std::vector<double_array> &w, int *Ns,double ds,double_array &k_vector,double_array &dxyz){

  int         i,j,l,s;
  double      Q;

  double_array 	qB4(Nx,Ny,Nz,(Ns[5]+1));
  double_array  qintB4(Nx,Ny,Nz);

  //+++++++++++++++++++++++++++++++++++++++Forward++++++++++++++++++++++++++++++++++++++++++++++++++
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(l=0;l<Nz;l++){
	qintB4(i,j,l)=1.0;
      }
    }
  }
  solveModDiffEqn_FFT_Forward(qB4,w[5],qintB4,ds,Ns[5],1,k_vector,dxyz);

  //++++++++++++++++++++++++++++++++++++++Single Chain Partition Function+++++++++++++++++++++++++++
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Q=0.0;
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(l=0;l<Nz;l++){
	Q+=(qB4(i,j,l,Ns[5])*dxyz(0)*dxyz(1)*dxyz(2));
      }
    }
  }
  // Normalizing with respect to the volume of the box
  Q/=((dxyz(0)*Nx)*(dxyz(1)*Ny)*(dxyz(2)*Nz));

  
  // Here we do the concentration calculation
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(l=0;l<Nz;l++){

	phi[5](i,j,l)=0.0;  //B4

	//B4
	for(s=0;s<(Ns[5]+1);s++){
	  if(s==0 || s==Ns[5]){
	    phi[5](i,j,l)+=0.5*qB4(i,j,l,s)*qB4(i,j,l,Ns[5]-s)*ds;
	  }else{
	    phi[5](i,j,l)+=qB4(i,j,l,s)*qB4(i,j,l,Ns[5]-s)*ds;
	  }
	}

	phi[5](i,j,l)*=(activity);

      }
    }
  }
  

  return Q;


};
