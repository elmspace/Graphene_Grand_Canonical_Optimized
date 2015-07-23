//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//                                Solving The Modified Diffusion Equation
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void solveModDiffEqn_FFT(double_array &q, double_array &w, double_array &qint, double ds, int Ns, int sign, double_array &k, double_array &dxyz){
  
  int            i,j,l,s,ss;  // some counters
  unsigned long  ijl; // This is used for the Fourier Transform
  double_array   wds(Nx,Ny,Nz);
  double_array   kds(Nx,Ny,Nz);

  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(l=0;l<Nz;l++){
	kds(i,j,l)=exp((-ds)*k(i,j,l));
	wds(i,j,l)=exp((-0.5)*ds*w(i,j,l));
      }
    }
  }
  
  if(sign==1){
   
    for(i=0;i<Nx;i++){
      for(j=0;j<Ny;j++){
	for(l=0;l<Nz;l++){
	  (q(i,j,l,0))=(qint(i,j,l));
	}
      }
    }
    
    for(s=0;s<Ns;s++){
      
      for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
	  for(l=0;l<Nz;l++){
	    ss=l+Nz*(j+Ny*i);
	    input_q[ss]=q(i,j,l,s)*wds(i,j,l);
	  }
	}
      }
      fftw_execute(forward_plan);


      for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
	  for(l=0;l<Nz;l++){
	    ss=l+Nz*(j+Ny*i);
	    transformed_q[ss]*=kds(i,j,l);
	  }
	}
      }
      fftw_execute(inverse_plan);
      
      for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
	  for(l=0;l<Nz;l++){
	    ss=l+Nz*(j+Ny*i);
	    q(i,j,l,s+1)=((final_q[ss]*wds(i,j,l))/(8.0*Nx*Ny*Nz));
	  }
	}
      }
    }

  }else{

    for(i=0;i<Nx;i++){
      for(j=0;j<Ny;j++){
	for(l=0;l<Nz;l++){
	  q(i,j,l,0)=qint(i,j,l);
	}
      }
    }

    for(s=0;s<(Ns);s++){
      
      for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
	  for(l=0;l<Nz;l++){
	    ss=l+Nz*(j+Ny*i);
	    input_q[ss]=q(i,j,l,s)*wds(i,j,l);
	  }
	}
      }
      fftw_execute(forward_plan);

      for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
	  for(l=0;l<Nz;l++){
	    ss=l+Nz*(j+Ny*i);
	    transformed_q[ss]*=kds(i,j,l);
	  }
	}
      }
      fftw_execute(inverse_plan);
      
      for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
	  for(l=0;l<Nz;l++){
	    ss=l+Nz*(j+Ny*i);
	    q(i,j,l,s+1)=((final_q[ss]*wds(i,j,l))/(8.0*Nx*Ny*Nz));
	  }
	}
      }
    }
  }


};
