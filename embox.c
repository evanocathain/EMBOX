/**********************************************
 * EMBOX: A Code to Solve Maxwell's Equations *
 * Aeuthot: Evan Keane                        *
 * Date: 09/02/2010                           *
 *                                            *
 * curl(B) = dE/dt + mu*J                     *
 * curl(E) = -dB/dt                           *
 * div(E) = rho/epsilon                       *
 * div(B) = 0                                 *
 *                                            *
 ******************************************** */ 

// header files
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<time.h>
#include"grid.h"
#include"particles.h"
#include"update.h"
#include"embox_funcs.h"
#include"initialise.h"

// constants
#define CSQUARED 8.98755179e16
// functions

// Main function
int main(int argc, char **argv)
{
  // set default values for VARS
  int NPARTICLES = 10, SIZE=10, NSTEPS=100000;

  int i,j,k,n,x_pos,y_pos,z_pos;
  int dump_fields=0,dump_posns=1;
  int mode=0;
  double dt=1e-12,dx=0.03,dy=0.03,dz=0.03;
  double B0=1.0e-1;
  double mu0=4.0*M_PI*10.0e-7,epsilon0=1.0/(mu0*CSQUARED),q=1.6e-19;

  //struct particles charges[NPARTICLES];
   //struct grid fields[SIZE][SIZE][SIZE];

  srand(time(NULL));
/*
  // Check command line arguments
  if (argc != 2 ) {
    fprintf(stderr,"Usage: embox <mode>\n\n");
    fprintf(stderr,"mode - 1=random, 2=on a circle\n");
    exit(-1);
  }
  mode=atoi(argv[1]);
*/

  // Parse command line arguments, modify values if found
  for (i=1; i<argc; i++){
    if ( strcmp(argv[i], "-mode") == 0 ){
        i++;
        mode = atoi(argv[i]);
    }
    else if ( strcmp(argv[i], "-size") == 0){
        i++;
        SIZE = atoi(argv[i]);
    }
    else if ( strcmp(argv[i], "-np") == 0){
        i++;
        NPARTICLES = atoi(argv[i]);
    }
    else if ( strcmp(argv[i], "-ns") == 0){
        i++;
        NSTEPS = atoi(argv[i]);
    }
  }


  // MALLOC the arrays of structs
  struct particles *charges = malloc(sizeof *charges * NPARTICLES);
  if (charges == NULL){
    printf("Error allocating memory for charges (Error 1)\n");
    exit(1);
  }

  struct grid ***fields = malloc(SIZE * sizeof *fields );
  if (fields == NULL ) {
    printf("Error allocating memory for fields (Error 1)\n");
    exit(1);
  }
  for (i=0; i < SIZE; i++){

    // try to malloc second dimension!
    fields[i] = malloc(SIZE * sizeof *fields[i]);
    if ( fields[i] == NULL ){
      printf("Error allocating memory for fields (Error 2)\n");
      exit(1);
    }

    for (j=0; j<SIZE; j++) {
      fields[i][j] = malloc(SIZE * sizeof *fields[i][j]);
      if ( fields[i][j] == NULL ) {
        printf("Error allocating memory for fields (Error 3)\n");
        exit(1);
      }
    }
  }


  /* Make a charge distribution */
  /* DEPENDING ON MODE SWITCH VALUE */
    if (mode==1) {
        initialise_box(charges, NPARTICLES, SIZE);
    }
    else if (mode==2) {
        initialise_sphere(charges, NPARTICLES, SIZE, dx);
    }
    else {
      fprintf(stderr,"No choice made for initial particle positions and velocities.\nUse \"-mode [1|2]\" to select simulation mode\nExiting.\n");
      exit(-2);
    }
  

  /* Initialise Fields */
  for(i=0; i<SIZE; i++){
    for(j=0; j<SIZE; j++){
      for(k=0; k<SIZE; k++){
	fields[i][j][k].E[0]=1.0e5;
	fields[i][j][k].E[1]=0.0;
	fields[i][j][k].E[2]=0.0;
	fields[i][j][k].B[0]=0.0;
	fields[i][j][k].B[1]=1.0e-2;
	fields[i][j][k].B[2]=0.0;
	fields[i][j][k].J[0]=0.0;
	fields[i][j][k].J[1]=0.0;
	fields[i][j][k].J[2]=0.0;
	fields[i][j][k].rho = 0.0;
	/*
        r=dx*sqrt(pow(i-offset,2)+pow(j-offset,2)+pow(k-offset,2));
	printf("r = %lf\n",r);
	fields[i][j][k].E[0]=0.0;
	fields[i][j][k].E[1]=0.0;
	fields[i][j][k].E[2]=0.0;
	fields[i][j][k].B[0]=B0*(pow(radius,3)/pow(r,5))*3.0*(i-offset)*dx*(k-offset)*dx;
	fields[i][j][k].B[1]=B0*(pow(radius,3)/pow(r,5))*3.0*(j-offset)*dx*(k-offset)*dx;
	fields[i][j][k].B[2]=B0*(pow(radius,3)/pow(r,5))*(3.0*pow((k-offset)*dx,2)-pow(r,2));
	printf("Bx= %.15lf By= %.15lf Bz= %.15lf\n",fields[i][j][k].B[0],fields[i][j][k].B[1],fields[i][j][k].B[2]);
	fields[i][j][k].J[0]=0.0;
	fields[i][j][k].J[1]=0.0;
	fields[i][j][k].J[2]=0.0;
	*/
	if (dump_fields==1){
	  printf("%d %d %d %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf\n",i,j,k,fields[i][j][k].E[0],fields[i][j][k].E[1],fields[i][j][k].E[2],fields[i][j][k].B[0],fields[i][j][k].B[1],fields[i][j][k].B[2],fields[i][j][k].J[0],fields[i][j][k].J[1],fields[i][j][k].J[2]);
	}
      }
    }
    if (dump_fields==1){
      printf("\n");
    }
  }      
  
  /* Time Evolution Loop                         *
   * 1. calculate rho & J                        *
   * 2. update E & B                             *
   * 3. calculate Lorentz force on each particle *
   * 4. move each particle                       */
  for(n=0;n<NSTEPS;n++){
    /* 1. Calculate rho & J */
    // need to zero rho and J everywhere before calculating each time
    // maybe do it as last thing in time evo. loop somehow
    for(i=0; i<NPARTICLES; i++){
      x_pos=(int)charges[i].x[0];
      y_pos=(int)charges[i].x[1];
      z_pos=(int)charges[i].x[2];
      fields[x_pos][y_pos][z_pos].rho += q;
      fields[x_pos][y_pos][z_pos].J[0]+=q*charges[i].u[0];
      fields[x_pos][y_pos][z_pos].J[1]+=q*charges[i].u[1];
      fields[x_pos][y_pos][z_pos].J[2]+=q*charges[i].u[2];
    }
    
    // use rho & J to calculate update E and B fields using the curl equations
    // loop over i,j,k inside the function
    updatefield(fields, SIZE, dx, dy, dz, dt, dump_fields);
    
    // calculate the corresponding Lorentz force on each particle
    updatecharges(charges, fields, NPARTICLES, dt, dump_posns);
    
    // Zero the rho and J values in "fields"
    resetfield_rho_j(fields, SIZE);
	
  } /* END OF TIME EVOLUTION LOOP */
  
  
  // FREE MEMORY JUST TO BE TIDY
  free_grid(fields, SIZE, SIZE);
  free(charges);
    
  return(0); /* THE END! */
}
/*
double pick_rand(void)
{
  double num;
  num=(double)rand();
  num=num/(RAND_MAX+0.0);
  return(num); // a number between 0 and 1
}

void free_grid(struct grid ***data, size_t xlen, size_t ylen)
{
    size_t i, j;
    for (i=0; i<xlen; i++){
        if (data[i] != NULL){
            for (j=0; j<ylen; j++){
                free(data[i][j]);
            }
            free(data[i]);
        }
    }
    free(data);
}
*/
