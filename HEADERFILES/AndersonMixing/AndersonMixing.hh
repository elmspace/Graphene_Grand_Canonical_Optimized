/*
  Here we have the Anderson Mixing algorithm, which can found at:
  J. Chem. Phys. Vol 120, No 1, Jan 2004
*/
using namespace std;
void AndersonMixing(std::vector<double_array> &w, std::vector<double_array> &newW, std::vector<double_array> &delW,double_array &delphi, double_array &dxyz){

  int h, i, j, k, s;
  int n, m, chain;
  int error;
  double **U_Matrix;
  double **U_inv_Matrix;
  double *V_vector;
  double *A_vector;
  double Matrix_array[(History-1)*(History-1)];


  U_Matrix = create_2d_double_array(History-1,History-1,"U_Matrix");
  U_inv_Matrix = create_2d_double_array(History-1,History-1,"U_inv_Matrix");
  A_vector = create_1d_double_array(History-1,"A_vector");
  V_vector = create_1d_double_array(History-1,"V_vector");


  // ________________________________________________________________________________
  if(Anderson==1){
    // Set everything to zero on the first iteration

    if(iter==0){
      for(chain=0;chain<ChainType;chain++){
	for(h=0;h<History;h++){
	  for(i=0;i<Nx;i++){
	    for(j=0;j<Ny;j++){
	      for(k=0;k<Nz;k++){
		DW[h][chain][i][j][k] = 0.0;
		W[h][chain][i][j][k] = 0.0;
	      }
	    }
	  }
	}
      }
    }

    for(chain=0;chain<ChainType;chain++){
      for(h=(History-1);h>0;h--){
	for(i=0;i<Nx;i++){
	  for(j=0;j<Ny;j++){
	    for(k=0;k<Nz;k++){
	      DW[h][chain][i][j][k] = DW[h-1][chain][i][j][k];
	      W[h][chain][i][j][k] = W[h-1][chain][i][j][k];
	    }
	  }
	}
      }
    }

    for(chain=0;chain<ChainType;chain++){
      for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
	  for(k=0;k<Nz;k++){
	    DW[0][chain][i][j][k] = newW[chain](i,j,k)-w[chain](i,j,k);
	    W[0][chain][i][j][k] = w[chain](i,j,k);
	  }
	}
      }
    }
    

    if(abs(deltaW)<0.01){
      global_index=1;
    }
    if((global_index==1) && (iter>100)){
 
    
    //if(iter>50){


      for(n=0;n<History-1;n++){
	for(m=0;m<History-1;m++){ 
	  U_Matrix[n][m] = 0.0;
	}
	V_vector[n] = 0.0;
      }
      
      
      for(n=0;n<(History-1);n++){
	for(m=0;m<(History-1);m++){
	  for(chain=0;chain<ChainType;chain++){

	    U_Matrix[n][m] += Inner_Prod_Func(DW[0][chain],DW[0][chain],dxyz);
	    U_Matrix[n][m] += Inner_Prod_Func(DW[n+1][chain],DW[m+1][chain],dxyz);
	    U_Matrix[n][m] -= Inner_Prod_Func(DW[0][chain],DW[m+1][chain],dxyz);
	    U_Matrix[n][m] -= Inner_Prod_Func(DW[n+1][chain],DW[0][chain],dxyz);
	    
	  }
	}
	for(chain=0;chain<ChainType;chain++){
	  
	  V_vector[n] += Inner_Prod_Func(DW[0][chain],DW[0][chain],dxyz);
	  V_vector[n] -= Inner_Prod_Func(DW[n+1][chain],DW[0][chain],dxyz);
    
	}
      }
      
      //___________________________ Matrix Inversion_____________________________________
      // Putting the Matrix into a 1D format
      for(n=0;n<History-1;n++){
	for(m=0;m<History-1;m++){
	  s=m+n*(History-1);
	  Matrix_array[s]=U_Matrix[n][m];
	}
      }
      // Calling the matrix inversion function (it uses Lapack)
      error = matrix_invert((History-1), Matrix_array);
      // Checking for errors
      if(error!=0){
	if(error>0){
	  std::cout<<"The matrix is singular!"<<std::endl;
	  exit(0);
	}else{
	  std::cout<<"The matrix has an illegal value!"<<std::endl;
	  exit(0);
	}
      }
      // Reverting back to 2D format for the matrix
      n=0;
      m=0;
      for(s=0; s<((History-1)*(History-1)); s++){
	U_inv_Matrix[n][m]=Matrix_array[s];
	if(s%(History-1)==0 && s!=0){
	  n++;
	  m=0;
	}
	m++;
      }
      //___________________________________________________________________________________
      
      
      // Calculating A_vector Eq. 14
      for(n=0;n<(History-1);n++){
	A_vector[n]=0.0;
	for(m=0;m<(History-1);m++){
	  A_vector[n]+=U_inv_Matrix[n][m]*V_vector[m];
	}
      }

 
      for(chain=0;chain<ChainType;chain++){
	for(i=0;i<Nx;i++){
	  for(j=0;j<Ny;j++){
	    for(k=0;k<Nz;k++){
	      w[chain](i,j,k) = 0.0;
	    }
	  }
	}
      }
      
      epsilon_delomega_anderson = 1.0 - pow(0.9,iter);
      // Calculating the new omega fields
      for(chain=0;chain<ChainType;chain++){
	for(i=0;i<Nx;i++){
	  for(j=0;j<Ny;j++){
	    for(k=0;k<Nz;k++){


	      for(h=0;h<(History-1);h++){
		w[chain](i,j,k) += epsilon_delomega_anderson*A_vector[h]*(DW[h+1][chain][i][j][k]-DW[0][chain][i][j][k]);
	      }
	      w[chain](i,j,k) += epsilon_delomega_anderson*DW[0][chain][i][j][k];
	      
	      for(h=0;h<History-1;h++){
		w[chain](i,j,k) += A_vector[h]*(W[h+1][chain][i][j][k]-W[0][chain][i][j][k]);
	      }
	      w[chain](i,j,k) += W[0][chain][i][j][k];
	      
	    }
	  }
	}
      }



    }else{

      // We use simple mixing to build some history, before turning on Anderson Mixing
      for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
	  for(k=0;k<Nz;k++){
	    for(chain=0;chain<ChainType;chain++){
	      w[chain](i,j,k)+=(epsilon_delomega*delW[chain](i,j,k)-epsilon_delphi*delphi(i,j,k));
	    }
	  }
	}
      }

    }

  }else{

    // If Anderson Mixing is turned off, then simple mixing is used
    for(i=0;i<Nx;i++){
      for(j=0;j<Ny;j++){
	for(k=0;k<Nz;k++){
	  for(chain=0;chain<ChainType;chain++){
	    w[chain](i,j,k)+=(epsilon_delomega*delW[chain](i,j,k)-epsilon_delphi*delphi(i,j,k));
	  }
	}
      }
    }

  }







  destroy_2d_double_array(U_Matrix);
  destroy_2d_double_array(U_inv_Matrix);
  destroy_1d_double_array(A_vector);
  destroy_1d_double_array(V_vector);

}
