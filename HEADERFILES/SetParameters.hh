void parameters(double_array &chi,double &ds,int *Ns,double_array &dxyz,double_array &chiMatrix){
  
  int Ds;
  double xAB,xAC,xBC, xAhom, xChom, xBhom;
  double a;
  double Lx,Ly,Lz;

  // Setting chemical potentials
  mu_homo=-30000.0;
  mu_copo=0.0;

  // 0 = Read From File    1 = Make Structure
  Iomega=0;
  
  // Minimize with respect to box size (yes=1, No=0)
  box_min=0;

  // Setting the generic chi parameters
  xAB=80.0;
  xAC=80.0;
  xBC=80.0;
  xAhom=80.0;
  xChom=80.0;
  xBhom=0.0;

  // Degree of polymerization
  Ns[0]=10; // A
  Ns[1]=10; // C
  Ns[3]=NB_middle; // B2 Center
  Ns[2]=(100-(Ns[0]+Ns[1]+Ns[3]))/2; // B1
  Ns[4]=Ns[2]; // B3
  Ns[5]=100; // B4 (Homopolymer)

  // Total length
  Ds=Ns[0]+Ns[1]+Ns[2]+Ns[3]+Ns[4];

  // Step size along chain
  ds=1.0/Ds;
  
  // Setting box size
  Lx=4.0;
  Ly=4.0;
  Lz=4.0;

  // Initial Guess for box dimension
  if(CAC==1){ Lx=2.64; Ly=4.36; Lz=2.5;}
  if(ZnSc==1){ Lx=4.; Ly=4.; Lz=4.;}
  if(AlphaBN==1){ Lx=2.4; Ly=4.2; Lz=4.8;}
  if(AlphaBNBilayer==1){ Lx=2.4; Ly=4.2; Lz=4.8;}

 
  // dx, dy, dz step size
  dxyz(0)=Lx/Nx;
  dxyz(1)=Ly/Ny;
  dxyz(2)=Lz/Nz;
 
  //Setting the chain fractions
  fA=(double)Ns[0]/(double)Ds;// fA
  fC=(double)Ns[1]/(double)Ds;// fC
  fB1=(double)Ns[2]/(double)Ds;// fB1
  fB2=(double)Ns[3]/(double)Ds;// fB2
  fB3=(double)Ns[4]/(double)Ds;// fB3
  kappa=(double)Ns[5]/(double)Ds;// kappa for homopolymer

  std::cout<<"fA= "<<fA<<std::endl;
  std::cout<<"fC= "<<fC<<std::endl;
  std::cout<<"fB1= "<<fB1<<std::endl;
  std::cout<<"fB2= "<<fB2<<std::endl;
  std::cout<<"fB3= "<<fB3<<std::endl;
  
  // Setting up the individual chi values
  chi(0)=xAC;  //xAC
  chi(1)=xAB;  //xAB1
  chi(2)=xAB;  //xAB2
  chi(3)=xAB;  //xAB3
  chi(4)=xAhom;  //xAB4
  chi(5)=xBC;  //xCB1
  chi(6)=xBC;  //xCB2
  chi(7)=xBC;  //xCB3
  chi(8)=xChom;  //xCB4
  chi(9)=0.001;  //xB1B2
  chi(10)=0.001; //xB1B3
  chi(11)=xBhom; //xB1B4
  chi(12)=0.001; //xB2B3
  chi(13)=xBhom; //xB2B4
  chi(14)=xBhom; //xB3B4
  // Setting up the chi matrix in this case 8 by 8
  //1
  chiMatrix(0,0)=0.0;       
  chiMatrix(0,1)=chi(0);    
  chiMatrix(0,2)=chi(1);     
  chiMatrix(0,3)=chi(2);     
  chiMatrix(0,4)=chi(3);     
  chiMatrix(0,5)=chi(4);    
  //2
  chiMatrix(1,0)=chi(0);       
  chiMatrix(1,1)=0.0;    
  chiMatrix(1,2)=chi(5);     
  chiMatrix(1,3)=chi(6);     
  chiMatrix(1,4)=chi(7);     
  chiMatrix(1,5)=chi(8);    
  //3
  chiMatrix(2,0)=chi(1);      
  chiMatrix(2,1)=chi(5);     
  chiMatrix(2,2)=0.0;     
  chiMatrix(2,3)=chi(9);     
  chiMatrix(2,4)=chi(10);     
  chiMatrix(2,5)=chi(11);    
  //4
  chiMatrix(3,0)=chi(2);     
  chiMatrix(3,1)=chi(6); 
  chiMatrix(3,2)=chi(9);   
  chiMatrix(3,3)=0.0;    
  chiMatrix(3,4)=chi(12);   
  chiMatrix(3,5)=chi(13);    
  //5
  chiMatrix(4,0)=chi(3);      
  chiMatrix(4,1)=chi(7);    
  chiMatrix(4,2)=chi(10);    
  chiMatrix(4,3)=chi(12);    
  chiMatrix(4,4)=0.0;    
  chiMatrix(4,5)=chi(14);    
  //6
  chiMatrix(5,0)=chi(4);       
  chiMatrix(5,1)=chi(8);     
  chiMatrix(5,2)=chi(11);    
  chiMatrix(5,3)=chi(13);    
  chiMatrix(5,4)=chi(14);     
  chiMatrix(5,5)=0.0;    
 
};
