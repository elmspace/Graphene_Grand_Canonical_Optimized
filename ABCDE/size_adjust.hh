double size_adjust(std::vector<double_array> &w, std::vector<double_array> &phi, double_array &eta, int *Ns, double ds, double_array &k_vector, double_array &chi, double_array &dxyz, double_array &chiMatrix){

  int     i,j,k,l,ii,jj,ll,kk;
  double  delx,dely,delz;
  double  box;
  double_array dxyz_temp(3);
  double  *box_x,*box_y,*box_z,*box_fE;
  std::vector<double_array> w_temp;

  w_temp.reserve(ChainType);
  for(int chain=0; chain!=ChainType; chain++) w_temp.push_back(double_array(Nx, Ny, Nz));

  box_x=create_1d_double_array(26,"box_x");
  box_y=create_1d_double_array(26,"box_y");
  box_z=create_1d_double_array(26,"box_z");
  box_fE=create_1d_double_array(26,"box_fE");

  delx=0.01/Nx;
  dely=0.01/Ny;
  delz=0.01/Nz;
  //std::cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  l=0;
  for(i=-1;i<2;i++){
    for(j=-1;j<2;j++){
      for(k=-1;k<2;k++){
	
	if((i==0)&&(j==0)&&(k==0)){
	  
	}else{
	  
	  box_x[l]=i*delx;
	  box_y[l]=j*dely;
	  box_z[l]=k*delz;
	  
	  dxyz_temp=dxyz;
	  
	  for(ii=0;ii<Nx;ii++){
	    for(jj=0;jj<Ny;jj++){
	      for(kk=0;kk<Nz;kk++){
		for(ll=0;ll<ChainType;ll++){
		  w_temp[ll](ii,jj,kk)=w[ll](ii,jj,kk);
		}}}}
	  
	  dxyz_temp(0)+=box_x[l];
	  dxyz_temp(1)+=box_y[l];
	  dxyz_temp(2)+=box_z[l];
	  
	  box_fE[l]=FreeEnergy_Box_Edition(w_temp,phi,eta,Ns,ds,k_vector,chi,dxyz_temp,chiMatrix);
	  
	  //std::cout<<l<<" "<<dxyz_temp[0]<<" "<<dxyz_temp[1]<<" "<<dxyz_temp[2]<<" "<<box_fE[l]<<std::endl;
	  
	  l=l+1;
	  
	}
	
      }
    }
  }
  //std::cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  
  for(i=0;i<25;i++){
    ii=0;
    for(j=i+1;j<26;j++){
      if(box_fE[i]<box_fE[j]){
	ii=ii+1;
      }
    }
    if(ii==(26-(i+1))){
      break;
    }
  }

  
  dxyz(0)=dxyz(0)+box_x[i];
  dxyz(1)=dxyz(1)+box_y[i];
  dxyz(2)=dxyz(2)+box_z[i];

  //std::cout<<i<<std::endl;
  
  
  box=1;
  return box;
  
  destroy_1d_double_array(box_x);
  destroy_1d_double_array(box_y);
  destroy_1d_double_array(box_z);
  destroy_1d_double_array(box_fE);

  
};

