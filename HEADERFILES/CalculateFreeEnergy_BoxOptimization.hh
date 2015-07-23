double FreeEnergy_Box_Edition(std::vector<double_array> &w_temp, std::vector<double_array> &phi, double_array &eta, int *Ns, double ds, double_array &k_vector, double_array &chi, double_array &dxyz_temp, double_array &chiMatrix){

  
  double  currentfE; 
  double  QMultiBlock,QHomo; 
  double  fES; 

  WaveVectors(k_vector,dxyz_temp);
  
  currentfE=0.0;
  fES=0.0;
  
  QMultiBlock=ConcMultiBlock(phi,w_temp,Ns,ds,k_vector,dxyz_temp);
  QHomo=ConcHomo(phi,w_temp,Ns,ds,k_vector,dxyz_temp);
    
  fES=QMultiBlock+activity*QHomo;
  
  currentfE=-fES;
  return currentfE;
  
  
};
