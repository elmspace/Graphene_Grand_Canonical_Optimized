/*
  Here we are calculating the inner product between two function, based on the definiotn from:
  J. Chem. Phys. Vol 120, No 1, Jan 2004 (Eq. 10 and Eq. 12, 13)
*/
double Inner_Prod_Func(double ***f1, double ***f2, double_array &dxyz){

  int i, j, k;
  double result = 0.0;

  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(k=0;k<Nz;k++){
	result += f1[i][j][k]*f2[i][j][k]*dxyz(0)*dxyz(1)*dxyz(2);
      }
    }
  }
  result/= (dxyz(0)*dxyz(1)*dxyz(2)*Nx*Ny*Nz);

  return result;
  
};
