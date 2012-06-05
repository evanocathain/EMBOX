#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#include"particles.h"
#include"grid.h"
#include"initialise.h"
#include"embox_funcs.h"

void initialise_distn_box(struct particles *charges,
			  int nparticles, 
			  int size,
			  double dx,
			  double dy,
			  double dz,
			  FILE *charges_fp){

    int i; // loop variable

    for (i=0; i<nparticles; i++){
        // initialise all particles
        charges[i].x[0]=(size)*dx*pick_rand();
        charges[i].x[1]=(size)*dy*pick_rand();
        charges[i].x[2]=(size)*dz*pick_rand();
        charges[i].u[0]=0;
        charges[i].u[1]=0;
        charges[i].u[2]=0;
        charges[i].q=1.0;
    }
}

void initialise_distn_sphere(struct particles *charges,
			     int nparticles,
			     int size,
			     double dx,
			     FILE *charges_fp){

    int i; // loop variable
    double offset = size*dx*0.5;
    double radius=dx,phi,theta;

    for (i=0; i<nparticles; i++){
        phi=(2*M_PI*pick_rand());   // between 0 and 2*M_PI
        theta=(M_PI*pick_rand());   // between 0 and +M_PI
        charges[i].x[0]=offset+radius*cos(phi)*sin(theta);
        charges[i].x[1]=offset+radius*sin(phi)*sin(theta);
        charges[i].x[2]=offset+radius*cos(theta);
        charges[i].u[0]=0.0;
        charges[i].u[1]=0.0;
        charges[i].u[2]=0.0;
	//        charges[i].q=1.0;
	if (pick_rand() >= 0.5){
	  charges[i].q=1.0;
	}else{
	  charges[i].q=-1.0;
	}
	fprintf(charges_fp,"%f\n",charges[i].q);
    }
}
