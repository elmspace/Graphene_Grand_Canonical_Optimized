int matrix_invert(int N, double *matrix){
  int error = 0;
  int *pivot = (int *)malloc(N*sizeof(int));
  double *workspace = (double *)malloc(N*sizeof(double));
  dgetrf_(&N, &N, matrix, &N, pivot, &error);

  if(error != 0){
    free(pivot);
    free(workspace);
    return error;
  }

  dgetri_(&N, matrix, &N, pivot, workspace, &N, &error);

  if(error!=0){
    free(pivot);
    free(workspace);
    return error;
  }

  free(pivot);
  free(workspace);
  return error;

}
