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
                    int size){

    int i; // loop variable

    for (i=0; i<nparticles; i++){
        // initialise all particles
        charges[i].x[0]=(size-1)*pick_rand();
        charges[i].x[1]=(size-1)*pick_rand();
        charges[i].x[2]=(size-1)*pick_rand();
        charges[i].u[0]=0;
        charges[i].u[1]=0;
        charges[i].u[2]=0;
        charges[i].q=1.0;
    }
}

void initialise_distn_sphere(struct particles *charges,
                       int nparticles,
                       int size,
                       double dx){

    int i; // loop variable
    int offset = size*0.5; // should offset be a double?
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
        charges[i].q=1.0;
    }
}
