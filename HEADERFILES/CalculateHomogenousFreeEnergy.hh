double homogenousfE(double_array &chiMatrix, double_array &chi){

  int i, j, iter = 0, Iter = 100000;
  double wA_ave, wC_ave, wB1_ave, wB2_ave, wB3_ave, wB4_ave;
  double pA_ave, pC_ave, pB1_ave, pB2_ave, pB3_ave, pB4_ave;
  double dwA_ave, dwC_ave, dwB1_ave, dwB2_ave, dwB3_ave, dwB4_ave;
  double eta_ave;
  double dW, dW_cutoff = 1.0e-3;
  double dP, dP_cutoff = 1.0e-3;
  double sigma = 0.005;
  double  *w_ave, *p_ave;
  double fE_homo;

  w_ave=create_1d_double_array(ChainType,"w_ave");
  p_ave=create_1d_double_array(ChainType,"p_ave");
  
  
  dwA_ave=0.0;
  dwC_ave=0.0;
  dwB1_ave=0.0;
  dwB2_ave=0.0;
  dwB3_ave=0.0;
  dwB4_ave=0.0;

  eta_ave=0.0;

  pA_ave=0.01;  
  pC_ave=pA_ave;
  pB1_ave=pA_ave;
  pB2_ave=pA_ave;
  pB3_ave=pA_ave;
  pB4_ave=1.0-(pA_ave+pC_ave+pB1_ave+pB2_ave+pB3_ave);


  wA_ave = chi(0)*pC_ave + chi(1)*pB1_ave + chi(2)*pB2_ave + chi(3)*pB3_ave + chi(4)*pB4_ave + eta_ave;
  wC_ave = chi(0)*pA_ave + chi(5)*pB1_ave + chi(6)*pB2_ave + chi(7)*pB3_ave + chi(8)*pB4_ave + eta_ave;
  wB1_ave = chi(1)*pA_ave + chi(5)*pC_ave + chi(9)*pB2_ave + chi(10)*pB3_ave + chi(11)*pB4_ave + eta_ave;
  wB2_ave = chi(2)*pA_ave + chi(6)*pC_ave + chi(9)*pB1_ave + chi(12)*pB3_ave + chi(13)*pB4_ave + eta_ave;
  wB3_ave = chi(3)*pA_ave + chi(7)*pC_ave + chi(10)*pB1_ave + chi(12)*pB2_ave + chi(14)*pB4_ave + eta_ave;
  wB4_ave = chi(4)*pA_ave + chi(8)*pC_ave + chi(11)*pB1_ave + chi(13)*pB2_ave + chi(14)*pB3_ave + eta_ave;

  do{

    eta_ave = eta_ave - 0.05*(1.0-(pA_ave+pC_ave+pB1_ave+pB2_ave+pB3_ave+pB4_ave));
    
    pA_ave = exp(mu_copo - wA_ave*fA - wC_ave*fC - wB1_ave*fB1 - wB2_ave*fB2 - wB3_ave*fB3) * fA;
    pC_ave = exp(mu_copo - wA_ave*fA - wC_ave*fC - wB1_ave*fB1 - wB2_ave*fB2 - wB3_ave*fB3) * fC;
    pB1_ave = exp(mu_copo - wA_ave*fA - wC_ave*fC - wB1_ave*fB1 - wB2_ave*fB2 - wB3_ave*fB3) * fB1;
    pB2_ave = exp(mu_copo - wA_ave*fA - wC_ave*fC - wB1_ave*fB1 - wB2_ave*fB2 - wB3_ave*fB3) * fB2;  
    pB3_ave = exp(mu_copo - wA_ave*fA - wC_ave*fC - wB1_ave*fB1 - wB2_ave*fB2 - wB3_ave*fB3) * fB3;

    pB4_ave = exp(mu_homo*kappa - wB4_ave*kappa);

    //std::cout<<pA_ave<<" "<<pB4_ave<<std::endl;

    dwA_ave = (chi(0)*pC_ave + chi(1)*pB1_ave + chi(2)*pB2_ave + chi(3)*pB3_ave + chi(4)*pB4_ave+ eta_ave) - wA_ave;
    dwC_ave = (chi(0)*pA_ave + chi(5)*pB1_ave + chi(6)*pB2_ave + chi(7)*pB3_ave + chi(8)*pB4_ave+ eta_ave) - wC_ave;
    dwB1_ave = (chi(1)*pA_ave + chi(5)*pC_ave + chi(9)*pB2_ave + chi(10)*pB3_ave + chi(11)*pB4_ave+ eta_ave) - wB1_ave;
    dwB2_ave = (chi(2)*pA_ave + chi(6)*pC_ave + chi(9)*pB1_ave + chi(12)*pB3_ave + chi(13)*pB4_ave+ eta_ave) - wB2_ave;
    dwB3_ave = (chi(3)*pA_ave + chi(7)*pC_ave + chi(10)*pB1_ave + chi(12)*pB2_ave + chi(14)*pB4_ave+ eta_ave) - wB3_ave;
    dwB4_ave = (chi(4)*pA_ave + chi(8)*pC_ave + chi(11)*pB1_ave + chi(13)*pB2_ave + chi(14)*pB3_ave+ eta_ave) - wB4_ave;
    
    dP = 1.0 - (pA_ave+pC_ave+pB1_ave+pB2_ave+pB3_ave+pB4_ave);
    dW = (dwA_ave + dwC_ave + dwB1_ave + dwB2_ave + dwB3_ave + dwB4_ave) / 6.0;

    //std::cout<<dW<<" "<<dP<<std::endl;
    
    wA_ave += sigma*dwA_ave-sigma*dP;
    wC_ave += sigma*dwC_ave-sigma*dP;
    wB1_ave += sigma*dwB1_ave-sigma*dP;
    wB2_ave += sigma*dwB2_ave-sigma*dP;
    wB3_ave += sigma*dwB3_ave-sigma*dP;
    wB4_ave += sigma*dwB4_ave-sigma*dP;

    iter++;
    
  }while(iter<Iter);
  
  if(abs(dW)>0.0001){
    std::cout<<"The homogenous free energy did not converege!"<<std::endl;
    exit(-1);
  }
  
  w_ave[0] = wA_ave;
  w_ave[1] = wC_ave;
  w_ave[2] = wB1_ave;
  w_ave[3] = wB2_ave;
  w_ave[4] = wB3_ave;
  w_ave[5] = wB4_ave;

  p_ave[0] = pA_ave;
  p_ave[1] = pC_ave;
  p_ave[2] = pB1_ave;
  p_ave[3] = pB2_ave;
  p_ave[4] = pB3_ave;
  p_ave[5] = pB4_ave;

  fE_homo = 0.0;
  for(i=0;i<ChainType;i++){
    for(j=i;j<ChainType;j++){
      
      fE_homo+=p_ave[i]*p_ave[j]*chiMatrix(i,j);
      
    }
    fE_homo-=p_ave[i]*w_ave[i];
  }

  fE_homo-=exp(mu_copo - wA_ave*fA - wC_ave*fC - wB1_ave*fB1 - wB2_ave*fB2 - wB3_ave*fB3);
  fE_homo-=exp(mu_homo*kappa - wB4_ave*kappa)/kappa;


  Phi_Copo_Dis = pA_ave+pC_ave+pB1_ave+pB2_ave+pB3_ave;
  Phi_Homo_Dis = pB4_ave;

  
  destroy_1d_double_array(w_ave);
  destroy_1d_double_array(p_ave);
  
  return fE_homo;

  

};
