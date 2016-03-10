#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct Hilbert{
  int npoints;
  double* rotation_table;
  double* sense_table;
  double*** quad_table;
  double**  quad_pos;
  double**  corner;
  double* x_pos;
  double* y_pos;
  int levgrid;
  int levscan;
  int side;
  double* key_list;

} hilbert;

typedef struct Morton{
  
  double* x_pos;
  double* y_pos;
  double** quad;
  double** corner;
  int levgrid;
  int levscan;
  int size;
  int npoints;
  int* key_list;

} curve;

void bit_set(int* a, int i){

  *a |= 1 << i;
  
}

int bit_check(int* a, int i){

  return *a >> i & 1;

}

void bit_shift(int* a, int i){
  (*a) <<= i;
}

void bit_clean(int* a, int i){

  (*a) &= ~(1 << i);
  
}

int extract_bits(int a, int start, int count){

  return (a >> start) | ((1 << count) - 1);
}

void morton_init(struct Morton* m, int scan, int size){

  if(scan <= 0){
    printf("\n\twrong levscan value\n");
    return;
  }

  int dim = (int) pow(4, scan);
  double* tmp = (double*) malloc(8 * sizeof(double));

  (*m).npoints = dim;
  (*m).x_pos = (double*) malloc(dim * sizeof(double));
  (*m).y_pos = (double*) malloc(dim * sizeof(double));
  (*m).key_list = (int*) malloc(dim * sizeof(int));

  (*m).levgrid = 0;
  (*m).levscan = scan;
  (*m).size = size;

  (*m).quad = (double**) malloc(2 * sizeof(double*));
  (*m).corner = (double**) malloc(2 * sizeof(double*));

  (*m).quad[0] = &tmp[0];
  (*m).quad[1] = &tmp[4];

  (*m).quad[0][0] = 0.5;
  (*m).quad[1][0] = 0.5;

  (*m).quad[0][1] = 1.5;
  (*m).quad[1][1] = 0.5;

  (*m).quad[0][2] = 0.5;
  (*m).quad[1][2] = 1.5;

  (*m).quad[0][3] = 1.5;
  (*m).quad[1][3] = 1.5;

  tmp = (double*) malloc(8 * sizeof(double));

  (*m).corner[0] = &tmp[0];
  (*m).corner[1] = &tmp[4];

  (*m).corner[0][0] = 0.;
  (*m).corner[1][0] = 0.;

  (*m).corner[0][1] = 1.;
  (*m).corner[1][1] = 0.;

  (*m).corner[0][2] = 0.;
  (*m).corner[1][2] = 1.;

  (*m).corner[0][3] = 1.;
  (*m).corner[1][3] = 1.;

}

void hilbert_init(struct Hilbert* h, int side, int levscan){

  int i, j, k, npt;
  double* tmp = (double*) malloc(8 * sizeof(double));
  double** tr = (double**)malloc(8 * sizeof(double*));
  
  npt = (int)pow(4, levscan);
  (*h).side = side;
  (*h).levscan = levscan;
  (*h).levgrid = 0;
  (*h).npoints = npt;
  (*h).rotation_table = (double*)malloc(4 * sizeof(double));
  (*h).sense_table = (double*)malloc(4 * sizeof(double));
  (*h).x_pos = (double*) malloc(npt * sizeof(double));
  (*h).y_pos = (double*) malloc(npt * sizeof(double));
  
  (*h).rotation_table[0] = 3.;
  (*h).rotation_table[1] = 0.;
  (*h).rotation_table[2] = 0.;
  (*h).rotation_table[3] = 1.;

  (*h).sense_table[0] = -1.;
  (*h).sense_table[1] = 1.;
  (*h).sense_table[2] = 1.;
  (*h).sense_table[3] = -1.;
  
  (*h).quad_pos = (double**) malloc(2 * sizeof(double*));
  (*h).corner = (double**) malloc(2 * sizeof(double*));

  (*h).quad_pos[0] = &tmp[0];
  (*h).quad_pos[1] = &tmp[4];

  (*h).quad_pos[0][0] = 0.5;
  (*h).quad_pos[1][0] = 0.5;

  (*h).quad_pos[0][1] = 0.5;
  (*h).quad_pos[1][1] = 1.5;

  (*h).quad_pos[0][2] = 1.5;
  (*h).quad_pos[1][2] = 1.5;

  (*h).quad_pos[0][3] = 1.5;
  (*h).quad_pos[1][3] = 0.5;

  tmp = (double*) malloc(8 * sizeof(double));

  (*h).corner[0] = &tmp[0];
  (*h).corner[1] = &tmp[4];

  (*h).corner[0][0] = 0.;
  (*h).corner[1][0] = 0.;

  (*h).corner[0][1] = 0.;
  (*h).corner[1][1] = 1.;

  (*h).corner[0][2] = 1.;
  (*h).corner[1][2] = 1.;

  (*h).corner[0][3] = 1.;
  (*h).corner[1][3] = 0.;

  tmp = (double*) malloc(16 * sizeof(double));

  (*h).quad_table = (double***) malloc(4 * sizeof(double***));

  for(j = 0; j < 8; j++)
    tr[j] = &tmp[2*j];
  
  for(j = 0; j < 4; j++){
    (*h).quad_table[j] = &tr[2*j];
    for(k = 0; k < 2; k++)
      (*h).quad_table[j][k] = &tmp[4*j + 2*k];
  }

  (*h).quad_table[0][0][0] = 0;
  (*h).quad_table[0][1][0] = 3;
  (*h).quad_table[0][0][1] = 1;
  (*h).quad_table[0][1][1] = 2;

  (*h).quad_table[1][0][0] = 1;
  (*h).quad_table[1][1][0] = 0;
  (*h).quad_table[1][0][1] = 2;
  (*h).quad_table[1][1][1] = 3;

  (*h).quad_table[2][0][0] = 2;
  (*h).quad_table[2][1][0] = 1;
  (*h).quad_table[2][0][1] = 3;
  (*h).quad_table[2][1][1] = 0;

  (*h).quad_table[3][0][0] = 3;
  (*h).quad_table[3][1][0] = 2;
  (*h).quad_table[3][0][1] = 0;
  (*h).quad_table[3][1][1] = 1;

  printf("\n");
  for(i = 0; i < 4; i++){
    for(k = 0; k < 2; k++)
      for(j = 0; j < 2; j++)
	printf("\tquad_table[%d][%d][%d] = %lg\n", i, j, k, (*h).quad_table[i][j][k]);
    printf("\n");
  }
}

void hilbert_dealloc(struct Hilbert* h){
  free((*h).x_pos);
  free((*h).y_pos);
  free((*h).rotation_table);
  free((*h).sense_table);
  free((*h).quad_table[0][0]);
  free((*h).quad_table[0]);
  free((*h).quad_table[1]);
  free((*h).quad_table[2]);
  free((*h).quad_table[3]);
  free((*h).quad_pos[0]);
  free((*h).corner[0]);

}

void morton_dealloc(struct Morton* m){
  free((*m).x_pos);
  free((*m).y_pos);
  free((*m).key_list);
  free((*m).quad[0]);
  free((*m).corner[0]);
  free((*m).quad);
  free((*m).corner);
}

void make_morton_grid(struct Morton* m, double* xgrid, double* ygrid, int n){

  int i, j, iad;

  if( (*m).levgrid + 1 >= (*m).levscan ){
    /* COPY RESULTS INTO X_GRID AND Y_GRID */

    printf("\n\tNpoints = %d, n = %d", (*m).npoints, n);
    
    for(i = 0; i < n; i++){
      (*m).x_pos[i] = xgrid[i];
      (*m).y_pos[i] = ygrid[i];
    }
    
    return;
  }
  double* xsub = (double*) malloc(4 * n * sizeof(double));
  double* ysub = (double*) malloc(4 * n * sizeof(double));
  
  (*m).levgrid++;

  for(i = 0; i < 4; i++){

    iad = i * n;

    for(j = 0; j < n; j++){
      xsub[iad + j] = (xgrid[j] + (*m).corner[0][i] * (*m).size ) * 0.5;
      ysub[iad + j] = (ygrid[j] + (*m).corner[1][i] * (*m).size ) * 0.5;
    }
  }
  make_morton_grid(m, xsub, ysub, 4*n);
  free(xsub);
  free(ysub);
}

void make_morton_key(struct Morton* m){

  int j, ix, iy, tmp, levmorton = 8, levkey = 0;

  double cj = pow(2, levmorton) / (*m).size;

  for(j = 0; j < (*m).npoints; j++){

    ix = (*m).x_pos[j] * cj;
    iy = (*m).y_pos[j] * cj;
    tmp = 0;

    while(levkey < levmorton){

      if(bit_check(&ix, levkey))
	bit_set(&tmp, 2*levkey);

      if(bit_check(&iy, levkey))
	bit_set(&tmp, 2*levkey+1);
      
      levkey++;
    }
    (*m).key_list[j] = tmp;
    levkey = 0;
  }

}

int main(){

  int i, j, test = 0, new_test;
  double *xgrid, *ygrid;
  FILE* out;
  struct Morton curve;
  struct Hilbert hilbert;
  
  xgrid = (double*) malloc(4 * sizeof(double));
  ygrid = (double*) malloc(4 * sizeof(double));

  morton_init(&curve, 4, 16);
  hilbert_init(&hilbert, 5, 16);
  printf("\n\tmorton curve lscan: %d\n\n", curve.levscan);
  for(i = 0; i < 2; i++){
    for(j = 0; j < 4; j++)
      printf("\t%lg", curve.quad[i][j]);
    printf("\n");
  }

  /* initializing xgrid and ygrid */
  for(j = 0; j < 4; j++){
    xgrid[j] = curve.quad[0][j] * 0.5 * curve.size;
    ygrid[j] = curve.quad[1][j] * 0.5 * curve.size;
  }

  make_morton_grid(&curve, xgrid, ygrid, 4);
  make_morton_key(&curve);

  printf("\n\n\tsaving output\n");
  out = fopen("curve.dat", "w");

  for(i = 0; i < curve.npoints; i++)
    fprintf(out, "%lg\t%lg\t%d\n", curve.x_pos[i], curve.y_pos[i], curve.key_list[i]);

  fclose(out);
  morton_dealloc(&curve);
  hilbert_dealloc(&hilbert);
  free(xgrid);
  free(ygrid);

  return 0;
}
