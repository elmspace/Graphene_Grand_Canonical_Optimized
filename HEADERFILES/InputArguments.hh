void inputArguments(int numb_of_args, char* arg_input[]){


  // Set all phase type to 0
  // Ionic crystals
  AlphaBN=0;
  AlphaBN_single=0;
  CAC=0;
  CAC_single=0;
  CsCl=0;
  CsCl_single=0;
  ZnSc=0;
  ZnSc_single=0;
  NaCl=0;
  NaCl_single=0;
  
  // Classical phases
  LAM=0;
  HEX=0;
  BCC=0;

  
  if(strcmp( arg_input[1], "ABN") == 0){
    AlphaBN=1;
    Phase_Type = "AlphaBN";
  }else if(strcmp( arg_input[1], "ABN_s") == 0){
    AlphaBN_single=1;
    Phase_Type = "AlphaBN_single";
  }else if(strcmp( arg_input[1], "CAC") == 0){
    CAC=1;
    Phase_Type = "CAC";
  }else if(strcmp( arg_input[1], "CAC_s") == 0){
    CAC_single=1;
    Phase_Type = "CAC_single";
  }else if(strcmp( arg_input[1], "ZnSc") == 0){
    ZnSc=1;
    Phase_Type = "ZnSc";
  }else if(strcmp( arg_input[1], "ZnSc_s") == 0){
    ZnSc_single=1;
    Phase_Type = "ZnSc_single";
  }else if(strcmp( arg_input[1], "NaCl") == 0){
    NaCl=1;
    Phase_Type = "NaCl";
  }else if(strcmp( arg_input[1], "NaCl_s") == 0){
    NaCl_single=1;
    Phase_Type = "NaCl_single";
  }else if(strcmp( arg_input[1], "CsCl") == 0){
    CsCl=1;
    Phase_Type = "CsCl";
  }else if(strcmp( arg_input[1], "CsCl_s") == 0){
    CsCl_single=1;
    Phase_Type = "CsCl_single";
  }
  else{
    std::cout<<"The phase you have chosen does not exists!"<<std::endl;
    exit(1);
  }


  NB_middle=atoi(arg_input[2]);

  

};
