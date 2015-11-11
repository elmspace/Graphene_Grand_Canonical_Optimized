/*
  Here we have the Anderson Mixing algorithm, which can found at:
  J. Chem. Phys. Vol 120, No 1, Jan 2004
*/
using namespace std;
void AndersonMixing(std::vector<double_array> &w, std::vector<double_array> &newW, std::vector<double_array> &delW, double_array &dxyz){

  int h, i, j, k;
  int n, m, chain;
  double **U_Matrix;
  double **U_inv_Matrix;
  double *V_vector;
  double *A_vector;

  
  U_Matrix = create_2d_double_array(History-1,History-1,"U_Matrix");
  U_inv_Matrix = create_2d_double_array(History-1,History-1,"U_inv_Matrix");
  A_vector = create_1d_double_array(History-1,"A_vector");
  V_vector = create_1d_double_array(History-1,"V_vector");

  
 
  // ________________________________________________________________________________

  
  // Set everything to zero on the first iteration
  if(global_index==0){
    for(h=0;h<History;h++){
      for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
	  for(k=0;k<Nz;k++){
	    DW_0[h][i][j][k]=0.0;
	    DW_1[h][i][j][k]=0.0;
	    DW_2[h][i][j][k]=0.0;
	    DW_3[h][i][j][k]=0.0;
	    DW_4[h][i][j][k]=0.0;
	    DW_5[h][i][j][k]=0.0;

	    W_0[h][i][j][k]=0.0;
	    W_1[h][i][j][k]=0.0;
	    W_2[h][i][j][k]=0.0;
	    W_3[h][i][j][k]=0.0;
	    W_4[h][i][j][k]=0.0;
	    W_5[h][i][j][k]=0.0;
	  }
	}
      }
    }
  }
 
  for(h=(History-1);h>0;h--){
    for(i=0;i<Nx;i++){
      for(j=0;j<Ny;j++){
	for(k=0;k<Nz;k++){
	  DW_0[h][i][j][k] = DW_0[h-1][i][j][k];
	  DW_1[h][i][j][k] = DW_1[h-1][i][j][k];
	  DW_2[h][i][j][k] = DW_2[h-1][i][j][k];
	  DW_3[h][i][j][k] = DW_3[h-1][i][j][k];
	  DW_4[h][i][j][k] = DW_4[h-1][i][j][k];
	  DW_5[h][i][j][k] = DW_5[h-1][i][j][k];

	  W_0[h][i][j][k] = W_0[h-1][i][j][k];
	  W_1[h][i][j][k] = W_1[h-1][i][j][k];
	  W_2[h][i][j][k] = W_2[h-1][i][j][k];
	  W_3[h][i][j][k] = W_3[h-1][i][j][k];
	  W_4[h][i][j][k] = W_4[h-1][i][j][k];
	  W_5[h][i][j][k] = W_5[h-1][i][j][k];
	}
      }
    }
  }
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(k=0;k<Nz;k++){
	DW_0[0][i][j][k] = newW[0](i,j,k)-w[0](i,j,k);
	DW_1[0][i][j][k] = newW[1](i,j,k)-w[1](i,j,k);
	DW_2[0][i][j][k] = newW[2](i,j,k)-w[2](i,j,k);
	DW_3[0][i][j][k] = newW[3](i,j,k)-w[3](i,j,k);
	DW_4[0][i][j][k] = newW[4](i,j,k)-w[4](i,j,k);
	DW_5[0][i][j][k] = newW[5](i,j,k)-w[5](i,j,k);

	W_0[0][i][j][k] = w[0](i,j,k);
	W_1[0][i][j][k] = w[1](i,j,k);
	W_2[0][i][j][k] = w[2](i,j,k);
	W_3[0][i][j][k] = w[3](i,j,k);
	W_4[0][i][j][k] = w[4](i,j,k);
	W_5[0][i][j][k] = w[5](i,j,k);
      }
    }
  }


  if(global_index>10){
    
    
    for(chain=0;chain<ChainType;chain++){
      
      // Calculating the U_matrix and V_vector Eq. 12 and 13
      for(n=0;n<History-1;n++){
	for(m=0;m<History-1;m++){
	  if(chain==0){U_Matrix[n][m] = Inner_Prod_Func(DW_0[0],DW_0[n+1],DW_0[0],DW_0[m+1],dxyz,1);}
	  if(chain==1){U_Matrix[n][m] = Inner_Prod_Func(DW_1[0],DW_1[n+1],DW_1[0],DW_1[m+1],dxyz,1);}
	  if(chain==2){U_Matrix[n][m] = Inner_Prod_Func(DW_2[0],DW_2[n+1],DW_2[0],DW_2[m+1],dxyz,1);}
	  if(chain==3){U_Matrix[n][m] = Inner_Prod_Func(DW_3[0],DW_3[n+1],DW_3[0],DW_3[m+1],dxyz,1);}
	  if(chain==4){U_Matrix[n][m] = Inner_Prod_Func(DW_4[0],DW_4[n+1],DW_4[0],DW_4[m+1],dxyz,1);}
	  if(chain==5){U_Matrix[n][m] = Inner_Prod_Func(DW_5[0],DW_5[n+1],DW_5[0],DW_5[m+1],dxyz,1);}
	}
	if(chain==0){V_vector[n] = Inner_Prod_Func(DW_0[0],DW_0[n+1],DW_0[0],DW_0[0],dxyz,0);}
	if(chain==1){V_vector[n] = Inner_Prod_Func(DW_1[0],DW_1[n+1],DW_1[0],DW_1[0],dxyz,0);}
	if(chain==2){V_vector[n] = Inner_Prod_Func(DW_2[0],DW_2[n+1],DW_2[0],DW_2[0],dxyz,0);}
	if(chain==3){V_vector[n] = Inner_Prod_Func(DW_3[0],DW_3[n+1],DW_3[0],DW_3[0],dxyz,0);}
	if(chain==4){V_vector[n] = Inner_Prod_Func(DW_4[0],DW_4[n+1],DW_4[0],DW_4[0],dxyz,0);}
	if(chain==5){V_vector[n] = Inner_Prod_Func(DW_5[0],DW_5[n+1],DW_5[0],DW_5[0],dxyz,0);}
      }
         
      //   *** Do the matrix inversion call here
      MatrixInversion(U_Matrix,History-1,U_inv_Matrix);
      
      
      // Calculating A_vector Eq. 14
      for(n=0;n<History-1;n++){
	A_vector[n]=0.0;
	for(m=0;m<History-1;m++){
	  A_vector[n]+=U_inv_Matrix[n][m]*V_vector[m];
	}
      }
      
      // Calculating the new omega fields
      for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
	  for(k=0;k<Nz;k++){
	    
	    w[chain](i,j,k) = 0.0;
	    for(h=0;h<History-1;h++){
	      if(chain==0){w[0](i,j,k) += A_vector[h]*(DW_0[h+1][i][j][k]-DW_0[0][i][j][k]);}
	      if(chain==1){w[1](i,j,k) += A_vector[h]*(DW_1[h+1][i][j][k]-DW_1[0][i][j][k]);}
	      if(chain==2){w[2](i,j,k) += A_vector[h]*(DW_2[h+1][i][j][k]-DW_2[0][i][j][k]);}
	      if(chain==3){w[3](i,j,k) += A_vector[h]*(DW_3[h+1][i][j][k]-DW_3[0][i][j][k]);}
	      if(chain==4){w[4](i,j,k) += A_vector[h]*(DW_4[h+1][i][j][k]-DW_4[0][i][j][k]);}
	      if(chain==5){w[5](i,j,k) += A_vector[h]*(DW_5[h+1][i][j][k]-DW_5[0][i][j][k]);}
	    }
	    if(chain==0){w[0](i,j,k) += DW_0[0][i][j][k];}
	    if(chain==1){w[1](i,j,k) += DW_1[0][i][j][k];}
	    if(chain==2){w[2](i,j,k) += DW_2[0][i][j][k];}
	    if(chain==3){w[3](i,j,k) += DW_3[0][i][j][k];}
	    if(chain==4){w[4](i,j,k) += DW_4[0][i][j][k];}
	    if(chain==5){w[5](i,j,k) += DW_5[0][i][j][k];}

	    w[chain](i,j,k) = epsilon_delomega*w[chain](i,j,k) ;


	    for(h=0;h<History-1;h++){
	      if(chain==0){w[0](i,j,k) += A_vector[h]*(W_0[h+1][i][j][k]-W_0[0][i][j][k]);}
	      if(chain==1){w[1](i,j,k) += A_vector[h]*(W_1[h+1][i][j][k]-W_1[0][i][j][k]);}
	      if(chain==2){w[2](i,j,k) += A_vector[h]*(W_2[h+1][i][j][k]-W_2[0][i][j][k]);}
	      if(chain==3){w[3](i,j,k) += A_vector[h]*(W_3[h+1][i][j][k]-W_3[0][i][j][k]);}
	      if(chain==4){w[4](i,j,k) += A_vector[h]*(W_4[h+1][i][j][k]-W_4[0][i][j][k]);}
	      if(chain==5){w[5](i,j,k) += A_vector[h]*(W_5[h+1][i][j][k]-W_5[0][i][j][k]);}
	    }
	    if(chain==0){w[0](i,j,k) += W_0[0][i][j][k];}
	    if(chain==1){w[1](i,j,k) += W_1[0][i][j][k];}
	    if(chain==2){w[2](i,j,k) += W_2[0][i][j][k];}
	    if(chain==3){w[3](i,j,k) += W_3[0][i][j][k];}
	    if(chain==4){w[4](i,j,k) += W_4[0][i][j][k];}
	    if(chain==5){w[5](i,j,k) += W_5[0][i][j][k];}

	    w[chain](i,j,k) = epsilon_delomega*w[chain](i,j,k) ;
	    
	  }
	}
      }
       
    }
    
  }else{

    
    for(i=0;i<Nx;i++){
      for(j=0;j<Ny;j++){
	for(k=0;k<Nz;k++){
	  for(chain=0;chain<ChainType;chain++){
	    w[chain](i,j,k)+=(epsilon_delomega*delW[chain](i,j,k));
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





