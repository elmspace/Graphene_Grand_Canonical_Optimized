void Incomp(double_array &eta, std::vector<double_array> &phi, double_array &delphi, double_array &chiMatrix, std::vector<double_array> &w){

  int     i,j,k,ii,jj;
  int     chain; 
  double  ptot;
  
  int error;
  int s;
  double Matrix_array[ChainType*ChainType];
  double chiMatrix_inv[ChainType][ChainType];
  
  double w_array[ChainType];
  double one_array[ChainType];
  double w_star_array[ChainType];
  double one_star_array[ChainType];
  double one_hat, w_hat;

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

  /*
  // **************************************************************************************************************
  // **************************************************************************************************************
  // **************************************************************************************************************


  // Finding the inverse of the chiMatrix &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  // Putting the Matrix into a 1D format
  for(ii=0;ii<ChainType;ii++){
    for(jj=0;jj<ChainType;jj++){
      s=jj+ii*(ChainType);
      Matrix_array[s]=chiMatrix(ii,jj);
    }
  }
  
  // Calling the matrix inversion function (it uses Lapack)
  error = matrix_invert(ChainType, Matrix_array);
  // Checking for errors
  if(error!=0){
    if(error>0){
      std::cout<<"The matrix is singular! (CalculateEta.hh)"<<std::endl;
      exit(0);
    }else{
      std::cout<<"The matrix has an illegal value! (CalculateEta.hh)"<<std::endl;
      exit(0);
    }
  }
  // Reverting back to 2D format for the matrix
  ii=0;
  jj=0;
  for(s=0; s<(ChainType*ChainType); s++){
    chiMatrix_inv[ii][jj]=Matrix_array[s];
    if(s%ChainType==0 && s!=0){
      ii++;
      jj=0;
    }
    jj++;
  }
  // &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  for(ii=0;ii<ChainType;ii++){
    one_array[ii] = 1.0;
  }
  
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(k=0;k<Nz;k++){

	one_hat = 0.0;
	w_hat = 0.0;
	
	for(ii=0;ii<ChainType;ii++){
	  w_array[ii] = w[ii](i,j,k);
	}
	
	for(ii=0;ii<ChainType;ii++){
	  for(jj=0;jj<ChainType;jj++){
	    w_star_array[ii] = chiMatrix_inv[ii][jj]*w_array[jj];
	    one_star_array[ii] = chiMatrix_inv[ii][jj]*one_array[jj];
	  }
	}

	for(ii=0;ii<ChainType;ii++){
	  one_hat += one_star_array[ii];
	  w_hat += w_star_array[ii];
	}

	eta(i,j,k) = (1.0/one_hat)*(w_hat-1.0);
	
      }
    }
  }

  // **************************************************************************************************************
  // **************************************************************************************************************
  // **************************************************************************************************************
  */

};
