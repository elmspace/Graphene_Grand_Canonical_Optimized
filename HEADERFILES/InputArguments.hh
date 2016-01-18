void inputArguments(int numb_of_args, char* arg_input[]){


  //Set all phase type to 0
  AlphaBN=0;
  AlphaBNBilayer=0;
  CAC=0;
  CsCl=0;
  ZnSc=0;
  Bilayer=0;
  LAM=0;
  HEX=0;
  BCC=0;
  
  if(strcmp( arg_input[1], "A") == 0){
    AlphaBN=1;
  }else if(strcmp( arg_input[1], "AB") == 0){
    AlphaBNBilayer=1;
  }else if(strcmp( arg_input[1], "Z") == 0){
    ZnSc=1;
  }else if(strcmp( arg_input[1], "C") == 0){
    CAC=1;
  }else if(strcmp( arg_input[1], "L") == 0){
    LAM=1;
  }else{
    std::cout<<"The phase you have chosen does not exists!"<<std::endl;
    exit(1);
  }


  NB_middle=atoi(arg_input[2]);

  

};
