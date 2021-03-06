void FreeEnergy(std::vector<double_array> &w, std::vector<double_array> &phi, double_array &eta, int *Ns, double ds, double_array &k_vector, double_array &chi, double_array &dxyz, double_array &chiMatrix){

  int     maxIter=250;
  int     i,j,k,chain,ii,jj;
  int     box_minimize = 1;
  double  currentfE, oldfE, deltafE, oldfE_iter;
  double  QMultiBlock,QHomo;
  double  fE_homo;
  double  fEW, fEchi, fES;
  double_array delphi(Nx,Ny,Nz);
  std::vector<double_array> delW;
  std::vector<double_array> newW;

  // Setting the acitivity
  activity=(1.0/kappa)*exp(kappa*mu_homo - mu_copo);

  delW.reserve(ChainType);
  newW.reserve(ChainType);

  for(unsigned n=0; n!=ChainType; n++) {
    delW.push_back(double_array(Nx,Ny,Nz));
    newW.push_back(double_array(Nx,Ny,Nz));
  }

  // Calculating the Homogenous Free Energy
  //fE_homo=homogenousfE(chiMatrix,chi);

  std::cout<<"Dis Copolymer Concentration:  "<<Phi_Copo_Dis<<"  Dis Homopolymer Concentration:   "<<Phi_Homo_Dis<<std::endl;

  oldfE=1.0e2;
  std::ofstream outputFile("./RESULTS/fE.dat");
  do{

    WaveVectors(k_vector,dxyz);

    currentfE=0.0;
    deltafE=0.0;
    iter=0;

    do{
      
      fEW=0.0;
      fEchi=0.0;
      fES=0.0;
      deltaW=0.0;

      QMultiBlock=ConcMultiBlock(phi,w,Ns,ds,k_vector,dxyz);
      QHomo=ConcHomo(phi,w,Ns,ds,k_vector,dxyz);
      Incomp(eta,phi,delphi,chiMatrix,w);
      Phi_Copo_Ord/=(Nx*Ny*Nz);
      Phi_Homo_Ord/=(Nx*Ny*Nz);


      for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
	  for(k=0;k<Nz;k++){

	    for(ii=0;ii<ChainType;ii++){
	      newW[ii](i,j,k)=0.0;
	      for(jj=0;jj<ChainType;jj++){

		newW[ii](i,j,k)+=chiMatrix(ii,jj)*phi[jj](i,j,k);
		fEchi+=phi[ii](i,j,k)*chiMatrix(ii,jj)*phi[jj](i,j,k)*dxyz(0)*dxyz(1)*dxyz(2);

	      }

	      newW[ii](i,j,k)+=eta(i,j,k);
	      fEW+=(newW[ii](i,j,k)*phi[ii](i,j,k)*dxyz(0)*dxyz(1)*dxyz(2));
	      
	    }

	  }
	}
      }

      /*
      std::cout<<newW[0](Nx/2,Ny/2,Nz/2)<<std::endl;
      CalculateAvrg(newW);
      std::cout<<newW[0](Nx/2,Ny/2,Nz/2)<<std::endl;
      */

      for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
	  for(k=0;k<Nz;k++){
	    for(ii=0;ii<ChainType;ii++){
	      delW[ii](i,j,k)=newW[ii](i,j,k)-w[ii](i,j,k);
	    }
	  }
	}
      }
	      
      deltaW = CalculateError(delW,w,dxyz);

      fEchi/=(2.0*((Nx*dxyz(0))*(Ny*dxyz(1))*(Nz*dxyz(2))));
      fEW/=(((Nx*dxyz(0))*(Ny*dxyz(1))*(Nz*dxyz(2))));
      fES=QMultiBlock+activity*QHomo;

      //currentfE=-fES-fEW+fEchi-fE_homo;
      currentfE=-fES-fEW+fEchi;

      deltafE=fabs(currentfE-oldfE_iter);

      oldfE_iter=currentfE;

      // Calculating the new omega fields
      //AndersonMixing(w,newW,delW,delphi,dxyz);
      SimpleMixing(w,newW,delW,delphi,dxyz);

      if(Test==1){
	std::cout<<"Iter="<<iter<<"   dfE="<<currentfE<<"   delW=" << deltaW<<"   pCopo="<<Phi_Copo_Ord<<"   pHom="<<Phi_Homo_Ord<<std::endl;
      }
      if(iter%100==0){SaveData(phi,w,dxyz,0);}

      global_iter = iter;
      iter++;
    }while((deltaW>precision) || (iter<maxIter));

    SaveData(phi,w,dxyz,0);

    outputFile <<currentfE<<" "<<fE_homo<<" "<<dxyz(0)*Nx<<" "<<dxyz(1)*Ny<<" "<<dxyz(2)*Nz<<std::endl;

    size_adjust(w,phi,eta,Ns,ds,k_vector,chi,dxyz,chiMatrix);

    if(oldfE<currentfE){
      box_minimize=0;
    }else{
      oldfE=currentfE;
    }
    if(box_min==0){ box_minimize=0;}

  }while(box_minimize==1);


  Free_Energy=currentfE+fE_homo;
  Free_Energy_Homo=fE_homo;
  Lx=dxyz(0)*Nx;
  Ly=dxyz(1)*Ny;
  Lz=dxyz(2)*Nz;


  outputFile <<"Done"<<std::endl;
  outputFile.close();


};
