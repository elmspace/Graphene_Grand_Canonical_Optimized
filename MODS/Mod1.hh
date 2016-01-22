/*

  This Mod will scan over the chemical potential, which starts from mu_homo<<0 will be phi_H = 0 to mu_homo>>0 will be 100% homopolymer

 */
void Mod1(std::vector<double_array> &w, std::vector<double_array> &phi, double_array &eta, int *Ns, double ds, double_array &k_vector, double_array &chi, double_array &dxyz, double_array &chiMatrix){

  double del_mu = 0.01;

  // Cleaning the .dat file
  std::ofstream outputFile37("./RESULTS/MOD1.dat");
  outputFile37 << std::endl;
  outputFile37.close();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // 1=on  0=off
  Test = 1;

  parameters(chi,ds,Ns,dxyz,chiMatrix);
  mu_homo=-15.0;
  mu_copo=0.0;
  activity=(1.0/kappa)*exp(kappa*(mu_homo - mu_copo));
  
  omega(w,chiMatrix);
  
  do{

    FreeEnergy(w,phi,eta,Ns,ds,k_vector,chi,dxyz,chiMatrix);

    std::cout<<kappa<<" "<<mu_homo<<" "<<Phi_Copo_Dis<<" "<<Phi_Homo_Dis<<" "<<Phi_Copo_Ord<<" "<<Phi_Homo_Ord<<" "<<Free_Energy<<" "<<Free_Energy_Homo<<" "<<Lx<<" "<<Ly<<" "<<Lz<<std::endl;
    
    std::ofstream outputFile37("./RESULTS/MOD1.dat" , ios::app);
    outputFile37 <<kappa<<" "<<mu_homo<<" "<<Phi_Copo_Dis<<" "<<Phi_Homo_Dis<<" "<<Phi_Copo_Ord<<" "<<Phi_Homo_Ord<<" "<<Free_Energy<<" "<<Free_Energy_Homo<<" "<<Lx<<" "<<Ly<<" "<<Lz<<std::endl;
    outputFile37.close();

    SaveData(phi,w,dxyz,1);
    
    if(Test==1){break;}
    
    mu_homo+=del_mu;
    mu_copo=0.0;
    activity=(1.0/kappa)*exp(kappa*(mu_homo - mu_copo));

  }while(Phi_Homo_Ord<0.9);
  
}


