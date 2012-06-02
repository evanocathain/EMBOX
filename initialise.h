#ifndef INITIALISE_H
#define INITIALISE_H

void initialise_distn_box(struct particles *charges,
                    int nparticles,
			  int size,
			  double dx,
			  double dy,
			  double dz);

void initialise_distn_sphere(struct particles *charges,
                       int nparticles,
                       int size,
                       double dx);

#endif
