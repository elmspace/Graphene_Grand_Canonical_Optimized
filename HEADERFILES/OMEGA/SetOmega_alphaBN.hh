void omega_alphaBN(double ****w){

  int i,j,k;
  double delta;

  delta = 4.8;
  
  k=0;
  std::cout<<k<<std::endl;
  //for(k=0;k<1;k++){
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      w[0][i][j][k]=-30.0*cos(i*2.0*Pi/Nx)*cos(j*2.0*Pi/Ny); //A
      w[1][i][j][k]=-30.0*cos(i*2.0*Pi/Nx)*cos((j*2.0*Pi/Ny)-delta); //C
    }
  }
    //}

 


  k=Nz/2;
  std::cout<<k<<std::endl;
  //for(k=9;k<10;k++){
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      w[0][i][j][k]=-30.0*cos(i*2.0*Pi/Nx)*cos((j*2.0*Pi/Ny)-delta); //A
      w[1][i][j][k]=-30.0*cos(i*2.0*Pi/Nx)*cos((j*2.0*Pi/Ny)); //C
    }
  }
    //}

  
  
  k=(Nz-1);
  std::cout<<k<<std::endl;
  //for(k=(Nz-1);k<Nz;k++){
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      w[0][i][j][k]=-30.0*cos(i*2.0*Pi/Nx)*cos(j*2.0*Pi/Ny); //A
      w[1][i][j][k]=-30.0*cos(i*2.0*Pi/Nx)*cos((j*2.0*Pi/Ny)-delta); //C
    }
  }
  //}

  k=(Nz/2)-(Nz/4);
  std::cout<<k<<std::endl;
  //for(k=(Nz-1);k<Nz;k++){
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      w[2][i][j][k]=-30.0; //B1
      w[4][i][j][k]=-30.0; //B3
    }
  }
  k=(Nz/2)+(Nz/4);
  std::cout<<k<<std::endl;
  //for(k=(Nz-1);k<Nz;k++){
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      w[2][i][j][k]=-30.0; //B1
      w[4][i][j][k]=-30.0; //B3
    }
  }


  //std::cin>>k;

};



