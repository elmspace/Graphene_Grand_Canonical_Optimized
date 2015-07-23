void omega_alphaBN_bilayer(double ****w){

  int i,j,k;
  double delta;

  delta = 4.2;
  
  for(k=((Nz/2)-3);k<((Nz/2)+3);k++){
    for(i=0;i<Nx;i++){
      for(j=0;j<Ny;j++){
	w[0][i][j][k]=-100.0*cos(i*2.0*Pi/Nx)*cos(j*2.0*Pi/Ny); //A
	w[1][i][j][k]=-100.0*cos(i*2.0*Pi/Nx)*cos((j*2.0*Pi/Ny)-delta); //C
      }
    }
  }

  for(k=((Nz/2)-3);k<((Nz/2)+3);k++){
    for(i=0;i<Nx;i++){
      for(j=0;j<Ny;j++){
	w[5][i][j][k]=30.0; //Homo
      }
    }
  }


};



