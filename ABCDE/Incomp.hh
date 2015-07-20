void Incomp(double_array &eta, std::vector<double_array> &phi, double_array &delphi){

  int     i,j,k;
  int     chain; 
  double  ptot;

  ptot=0.0;

  Phi_Copo_Ord=0.0;
  Phi_Homo_Ord=0.0;

  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(k=0;k<Nz;k++){
	
	ptot=0.0;
	delphi(i,j,k)=0.0;
    
	for(chain=0;chain<ChainType;chain++){
	  ptot+=phi[chain](i,j,k);
	}

	delphi(i,j,k)=1.0-ptot;
	eta(i,j,k)-=delphi(i,j,k);


	Phi_Copo_Ord+=phi[0](i,j,k)+phi[1](i,j,k)+phi[2](i,j,k)+phi[3](i,j,k)+phi[4](i,j,k);
	Phi_Homo_Ord+=phi[5](i,j,k);
	  
      }
    }
  }
 
  

};
