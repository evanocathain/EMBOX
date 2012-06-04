/**********************************************
 * EMBOX: A Code to Solve Maxwell's Equations *
 * Author: Evan Keane                         *
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
    int NPARTICLES = 10, SIZE=10, NSTEPS=10000;

    int i,j,k,n;
    int dump_fields=0,dump_posns=1;
    int mode=0;
    double dt=1e-11,dx=0.03,dy=0.03,dz=0.03;
    double r,x,y,z,radius=dx,offset;
    double B0=5.0e-1;
    double mu0=4.0*M_PI*10.0e-7,epsilon0=1.0/(mu0*CSQUARED);
    FILE *positions_fp,*charges_fp,*fields_fp;

    srand(time(NULL));

    // Parse command line arguments, modify values if found
    if (argc == 1){
      fprintf(stdout,"Usage: embox -mode [1|2] <options> \n\n");
      fprintf(stdout,"Options:\t -mode\t 1=uniform particle dist., 2=distributed on a sphere\n");
      fprintf(stdout,"\t\t -size\t dimenstion of box [default=10] \n");
      fprintf(stdout,"\t\t -np\t number of particles to simulate [default=10] \n");
      fprintf(stdout,"\t\t -ns\t number of integration steps [default=10000] \n\n");
    }
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
    if (mode == 0){ // mode hasn't been set. 
      fprintf(stderr,"No choice made for initial particle positions and velocities.\nUse \"-mode [1|2]\" to select simulation mode\nExiting.\n");
      exit(-2);
    }


    // MALLOC the arrays of structs  
    struct particles *charges = malloc(sizeof *charges * NPARTICLES);
    if (charges == NULL){
      fprintf(stderr,"Error allocating memory for charges (Error 1)\n");
      exit(1);
    }

    struct grid ***fields = malloc(SIZE * sizeof *fields );
    if (fields == NULL ) {
      fprintf(stderr,"Error allocating memory for fields (Error 1)\n");
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


    offset=SIZE*0.5;

  /* Make a charge distribution */
  /* DEPENDING ON MODE SWITCH VALUE */

    charges_fp = fopen("charges","w");

    if (mode==1) {
      initialise_distn_box(charges, NPARTICLES, SIZE, dx, dy, dz, charges_fp);
    }
    else if (mode==2) {
      initialise_distn_sphere(charges, NPARTICLES, SIZE, dx, charges_fp);
    }
    

    // Open file for outputting the initial field values
    fields_fp = fopen("fields","w");

  /* Initialise Fields */
    /* THIS SHOULDPROBABLY GO INTO INITIALISE.C ONCE WE'VE FIGURED OUT
     * A SENSIBLE FIELD CONFIGURATION*/
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
	///*
	x=(i-offset)*dx;
	y=(j-offset)*dy;
	z=(k-offset)*dz;
	if (i==offset && j==offset && k==offset){
	  r=0.01*dx;
	  x=y=z=0.0;
	  fields[i][j][k].E[0]=0.0;
	  fields[i][j][k].E[1]=0.0;
	  fields[i][j][k].E[2]=0.0;
	  fields[i][j][k].B[0]=0.0;
	  fields[i][j][k].B[1]=0.0;
	  fields[i][j][k].B[2]=0.0;
	}else {
	  r=sqrt(pow(x,2)+pow(y,2)+pow(z,2));
	  fields[i][j][k].E[0]=0.0;
	  fields[i][j][k].E[1]=0.0;
	  fields[i][j][k].E[2]=0.0;
	  fields[i][j][k].B[0]=B0*(pow(radius,3)/pow(r,5))*3.0*x*z;
	  fields[i][j][k].B[1]=B0*(pow(radius,3)/pow(r,5))*3.0*y*z;
	  fields[i][j][k].B[2]=B0*(pow(radius,3)/pow(r,5))*(3.0*pow(z,2)-pow(r,2));
	}
	fprintf(fields_fp,"Bx= %.15lf By= %.15lf Bz= %.15lf x= %lf y= %lf z= %lf\n",fields[i][j][k].B[0],fields[i][j][k].B[1],fields[i][j][k].B[2],x,y,z);
	fields[i][j][k].J[0]=0.0;
	fields[i][j][k].J[1]=0.0;
	fields[i][j][k].J[2]=0.0;
	//*/
	if (dump_fields==1){
	  printf("%d %d %d %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf\n",i,j,k,fields[i][j][k].E[0],fields[i][j][k].E[1],fields[i][j][k].E[2],fields[i][j][k].B[0],fields[i][j][k].B[1],fields[i][j][k].B[2],fields[i][j][k].J[0],fields[i][j][k].J[1],fields[i][j][k].J[2]);
	}
      }
    }
    if (dump_fields==1){
      printf("\n");
    }
  }      



  // Open file for outputting positions of particles
  positions_fp = fopen("positions","w");


  
  /* Time Evolution Loop                         *
   * 1. calculate rho & J                        *
   * 2. update E & B                             *
   * 3. calculate Lorentz force on each particle *
   * 4. move each particle                       */
    for(n=0;n<NSTEPS;n++){
        /* 1. Calculate rho & J */
        // need to zero rho and J everywhere before calculating each time
        // maybe do it as last thing in time evo. loop somehow
        update_field_current(charges, fields, NPARTICLES);
      
        // use rho & J to calculate update E and B fields using the curl equations
        // loop over i,j,k inside the function
        update_field_strength(fields, SIZE, dx, dy, dz, dt, dump_fields);
    
        // calculate the corresponding Lorentz force on each particle
        update_charge_posns(charges, fields, NPARTICLES, dt, dx, dy, dz, SIZE, dump_posns, positions_fp);
    
        // Zero the rho and J values in "fields"
        resetfield_rho_j(fields, SIZE);

    } /* END OF TIME EVOLUTION LOOP */
  
  
    // FREE MEMORY JUST TO BE TIDY
    free_grid(fields, SIZE, SIZE);
    free(charges);
    
    fclose(positions_fp);
    fclose(fields_fp);
    fclose(charges_fp);

    return(0); /* THE END! */
}

