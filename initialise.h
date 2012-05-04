#ifndef INITIALISE_H
#define INITIALISE_H

void initialise_box(struct particles *charges,
                    int nparticles,
                    int size);

void initialise_sphere(struct particles *charges,
                       int nparticles,
                       int size,
                       double dx);

#endif
