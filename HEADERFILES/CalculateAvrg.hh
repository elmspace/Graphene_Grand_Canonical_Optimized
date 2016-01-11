/*
  Here we calculate the average of the omega fields for a given iteration
*/
void CalculateAvrg(std::vector<double_array> &AvgW){

  int i, j, k, chain;
  double avrg = 0.0;

  for(chain=0;chain<ChainType;chain++){
    
    avrg = 0.0;
    for(i=0;i<Nx;i++){
      for(j=0;j<Ny;j++){
	for(k=0;k<Nz;k++){
	  avrg += AvgW[chain](i,j,k);
	}
      }
    }
    avrg/=(Nx*Ny*Nz);
    
    for(i=0;i<Nx;i++){
      for(j=0;j<Ny;j++){
	for(k=0;k<Nz;k++){
	  AvgW[chain](i,j,k) -= avrg;
	}
      }
    }
    
  }

  return;
  
};
