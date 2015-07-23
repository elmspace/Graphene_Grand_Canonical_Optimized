void omega_CAC(double ****w){

  int i,j,k;
  
  for(k=0;k<Nz;k++){
    for(i=0;i<Nx;i++){
      for(j=0;j<Ny;j++){
	w[0][i][j][k]=-30.0*cos(i*2.0*Pi/Nx)*cos(j*2.0*Pi/Ny)*sin(k*2.0*Pi/Nz); //A
	w[1][i][j][k]=-30.0*cos(i*2.0*Pi/Nx)*cos(j*2.0*Pi/Ny); //C
      }
    }
  }

 


};



