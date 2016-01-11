/*
  Here we calculate the error in omega fields for a given iteration
*/
double CalculateError(std::vector<double_array> &delW, std::vector<double_array> &w, double_array &dxyz){

  int i, j, k, chain;
  double dw_sum = 0.0;
  double w_sum = 0.0;

  
  
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(k=0;k<Nz;k++){

	for(chain=0;chain<ChainType;chain++){

	  dw_sum += pow(delW[chain](i,j,k),2)*dxyz(0)*dxyz(1)*dxyz(2);
	  w_sum += pow(w[chain](i,j,k),2)*dxyz(0)*dxyz(1)*dxyz(2);

	}
	
      }
    }
  }

  return sqrt(dw_sum/w_sum);
  
};
