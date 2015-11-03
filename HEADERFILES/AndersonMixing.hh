/*
  Here we have the Anderson Mixing algorithm, which can found at:
  J. Chem. Phys. Vol 120, No 1, Jan 2004
*/
using namespace std;
void AndersonMixing(std::vector<double_array> &w, std::vector<double_array> &newW, double_array &dxyz){

  int h, i, j, k;
  int n, m;
  double **U_Matrix;
  double **U_inv_Matrix;
  double *V_vector;
  double *A_vector;

  
  U_Matrix = create_2d_double_array(History,History,"U_Matrix");
  U_inv_Matrix = create_2d_double_array(History,History,"U_inv_Matrix");
  A_vector = create_1d_double_array(History,"A_vector");
  V_vector = create_1d_double_array(History,"V_vector");

  
  std::vector<double_array> DW;
  DW.reserve(History);
  for(int n=0; n<History; n++) 
  {
    DW.push_back(double_array(Nx,Ny,Nz));
  }

  // Set everything to zero on the first iteration
  if(global_index==0){
    for(h=0;h<History;h++){
      for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
	  for(k=0;k<Nz;k++){
	    DW[h](i,j,k)=0.0;
	  }
	}
      }
    }
  }
  
 
  for(h=(History-1);h>0;h--){
    for(i=0;i<Nx;i++){
      for(j=0;j<Ny;j++){
	for(k=0;k<Nz;k++){
	  DW[h](i,j,k)=DW[h-1](i,j,k);
	}
      }
    }
  }
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(k=0;k<Nz;k++){
	DW[0](i,j,k)=newW[0](i,j,k)-w[0](i,j,k);
      }
    }
  }


  
  // Calculating the U_matrix and V_vector Eq. 12 and 13
  for(n=0;n<History;n++){
    for(m=0;m<History;m++){
      U_Matrix[n][m] = Inner_Prod_Func(DW[0],DW[n],DW[0],DW[m],dxyz,1);
    }
    V_vector[n] = Inner_Prod_Func(DW[0],DW[n],DW[0],DW[0],dxyz,0);
  }


  for(n=0;n<History;n++){
    for(m=0;m<History;m++){
      cout<<n<<" "<<m<<" "<<U_inv_Matrix[n][m]<<endl;
    }
  }
  cin>>i;

  
  // Calculating A_vector Eq. 14
  for(n=0;n<History;n++){
    for(m=0;m<History;m++){
      A_vector[n]+=U_Matrix[n][m]*V_vector[m];
    }
  }
  
  
  //   *** Do the matrix inversion call here
  MatrixInversion(U_Matrix,History,U_inv_Matrix);
  
  // Calculating the new omega fields
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(k=0;k<Nz;k++){

	for(h=0;h<History;h++){
	  w[0](i,j,k)+=A_vector[h]*(DW[h](i,j,k)-DW[0](i,j,k));
	}
	w[0](i,j,k)+=newW[0](i,j,k); 
	
      }
    }
  }


  destroy_2d_double_array(U_Matrix);
  destroy_2d_double_array(U_inv_Matrix);
  destroy_1d_double_array(A_vector);
  destroy_1d_double_array(V_vector);
  
}
