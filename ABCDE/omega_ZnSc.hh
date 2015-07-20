void omega_ZnSc(double ****w){

  int i,j,k;
  double delta;

  delta = 4.8;
  

  k=0;
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){  
      w[0][i][j][k]=-50.0*cos(i*2.0*Pi/Nx)*cos(j*2.0*Pi/Ny)*cos(k*2.0*Pi/Nz); //A
      w[1][i][j][k]=50.0;
    }
  }
  
  k=Nz/2;
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      w[0][i][j][k]=50.0*cos(i*2.0*Pi/Nx)*cos(j*2.0*Pi/Ny)*cos(k*2.0*Pi/Nz); //A
      w[1][i][j][k]=50.0;
    }
  }

  k=Nz-1;
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      w[0][i][j][k]=-50.0*cos(i*2.0*Pi/Nx)*cos(j*2.0*Pi/Ny)*cos(k*2.0*Pi/Nz); //A
      w[1][i][j][k]=50.0;
    }
  } 




  
  w[1][Nx/4][Ny/4][Nz/4]=-300.0; //C
  w[1][3*Nx/4][Ny/4][3*Nz/4]=-300.0; //C
  w[1][3*Nx/4][3*Ny/4][Nz/4]=-300.0; //C
  w[1][Nx/4][3*Ny/4][3*Nz/4]=-300.0; //C
 

};



