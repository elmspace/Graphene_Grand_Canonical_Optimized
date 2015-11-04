/*
  Here we are calculating the inner product between two function, based on the definiotn from:
  J. Chem. Phys. Vol 120, No 1, Jan 2004 (Eq. 10 and Eq. 12, 13)
*/
double Inner_Prod_Func(double_array &f1, double_array &f2, double_array &f3, double_array &f4, double_array &dxyz, int cond){

  int i, j, k;
  double result;
  
  
  std::vector<double_array> g;
  g.reserve(2);
  for(int n=0; n<2; n++) 
  {
    g.push_back(double_array(Nx,Ny,Nz));
  }


  for(i=0; i<Nx; i++){
    for(j=0; j<Ny; j++){
      for(k=0; k<Nz; k++){
	g[0](i,j,k)=f1(i,j,k)-f2(i,j,k);
	if(cond==1){g[1](i,j,k)=f3(i,j,k)-f4(i,j,k);}
	else if(cond==0){g[1](i,j,k)=f3(i,j,k)-0.0;}
	else{std::cout<<"Something wrong in Inner_Prod_Func"<<std::endl; exit(1);}
      }
    }
  }

  result=0.0;

  for(i=0; i<Nx; i++){
    for(j=0; j<Ny; j++){
      for(k=0; k<Nz; k++){
	result+=g[0](i,j,k)*g[1](i,j,k)*dxyz(0)*dxyz(1)*dxyz(2);
      }
    }
  }
  
  // Setting to zero, if number is too small
  if(abs(result)<1.0e-20){
    result = 0.0;
  }


  return result;
  
};
