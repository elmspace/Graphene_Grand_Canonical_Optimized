/*

  This Mod will scan over the chemical potential, which starts from mu_homo<<0 will be phi_H = 0 to mu_homo>>0 will be 100% homopolymer

 */
void Mod1(std::vector<double_array> &w, std::vector<double_array> &phi, double_array &eta, int *Ns, double ds, double_array &k_vector, double_array &chi, double_array &dxyz, double_array &chiMatrix){

  double del_mu = 0.025;

  // Tag1 and Tag2 are used for naming the output files
  Tag2 = (long long)(NB_middle);
  
  // Cleaning specific output file
  std::string OutPut="./RESULTS/MOD1_"+Phase_Type+"_"+std::to_string(Tag2)+".dat";
  std::ofstream outputFile57(OutPut);
  outputFile57 << std::endl;
  outputFile57.close();
  
  
  // Cleaning Generic output file
  std::ofstream outputFile37("./RESULTS/MOD1.dat");
  outputFile37 << std::endl;
  outputFile37.close();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // 1=on  0=off
  Test = 0;

  parameters(chi,ds,Ns,dxyz,chiMatrix);
  mu_homo=-15.00;
  mu_copo=0.0;
  activity=(1.0/kappa)*exp(kappa*(mu_homo - mu_copo));
  
  omega(w,chiMatrix);

  do{

    // Tag1 and Tag2 are used for naming the output files
    Tag1 = (long long)(mu_homo*100);
    Tag2 = (long long)(NB_middle);

    
    FreeEnergy(w,phi,eta,Ns,ds,k_vector,chi,dxyz,chiMatrix);

    // Output to screen
    std::cout<<kappa<<" "<<mu_homo<<" "<<Phi_Copo_Dis<<" "<<Phi_Homo_Dis<<" "<<Phi_Copo_Ord<<" "<<Phi_Homo_Ord<<" "<<Free_Energy<<" "<<Free_Energy_Homo<<" "<<Lx<<" "<<Ly<<" "<<Lz<<std::endl;

    // *********************************************
    // Writing to the output file with specific name
    std::string OutPut="./RESULTS/MOD1_"+Phase_Type+"_"+std::to_string(Tag2)+".dat";
    std::ofstream outputFile57(OutPut, ios::app);
    outputFile57 <<kappa<<" "<<mu_homo<<" "<<Phi_Copo_Dis<<" "<<Phi_Homo_Dis<<" "<<Phi_Copo_Ord<<" "<<Phi_Homo_Ord<<" "<<Free_Energy<<" "<<Free_Energy_Homo<<" "<<Lx<<" "<<Ly<<" "<<Lz<<std::endl;
    outputFile57.close();
    // *********************************************

    // *********************************************
    // Writing data to a generic file
    std::ofstream outputFile37("./RESULTS/MOD1.dat" , ios::app);
    outputFile37 <<kappa<<" "<<mu_homo<<" "<<Phi_Copo_Dis<<" "<<Phi_Homo_Dis<<" "<<Phi_Copo_Ord<<" "<<Phi_Homo_Ord<<" "<<Free_Energy<<" "<<Free_Energy_Homo<<" "<<Lx<<" "<<Ly<<" "<<Lz<<std::endl;
    outputFile37.close();
    // *********************************************
    
    SaveData(phi,w,dxyz,1);
    
    if(Test==1){break;}
    
    mu_homo+=del_mu;
    mu_copo=0.0;
    activity=(1.0/kappa)*exp(kappa*(mu_homo - mu_copo));
    
  }while(Phi_Homo_Ord<0.999);
  
}


