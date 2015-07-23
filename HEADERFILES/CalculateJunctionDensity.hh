void Junc_Density(double_array &qA,double_array &qC,double_array &qB1,double_array &qB2,double_array &qB3,double_array &qdagA,double_array qdagC,double_array &qdagB1,double_array &qdagB2,double_array &qdagB3,int *Ns,double ds, double Q){

  int i,j,l,k;
  double_array JB1C(Nx,Ny,Nz);
  double_array JCB2(Nx,Ny,Nz);
  double_array JB2A(Nx,Ny,Nz);
  double_array JAB3(Nx,Ny,Nz);
  
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(l=0;l<Nz;l++){
	JB1C(i,j,l)=qB1(i,j,l,Ns[2])*qdagB1(i,j,l,0)*ds;
	JCB2(i,j,l)=qC(i,j,l,Ns[1])*qdagC(i,j,l,0)*ds;
	JB2A(i,j,l)=qB2(i,j,l,Ns[3])*qdagB2(i,j,l,0)*ds;	
	JAB3(i,j,l)=qA(i,j,l,Ns[0])*qdagA(i,j,l,0)*ds;

	JB1C(i,j,l)*=(1.0/Q);
	JCB2(i,j,l)*=(1.0/Q);
	JB2A(i,j,l)*=(1.0/Q);
	JAB3(i,j,l)*=(1.0/Q);
      }
    }
  }

  

  //+++++++++++++++++++++++++++++ This output is setup for the matlab plotting +++++++++++++++++++  
  std::ofstream outputFile38("./ABCD_Junc.dat");
  for (i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(k=0;k<Nz;k++){//format A, C, B1+B2+B3, B4
	outputFile38<<JB1C(i,j,k)<<" "<<JCB2(i,j,k)<<" "<<JB2A(i,j,k)<<" "<<JAB3(i,j,k)<<std::endl;
      }}}
  outputFile38.close();
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  

};
