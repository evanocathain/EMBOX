#ifndef UPDATE_H
#define UPDATE_H

void updatefield(struct grid ***fields,
                int size, 
                double dx, 
                double dy, 
                double dz, 
                double dt,
                int dump);

void updatecharges(struct particles *charges, 
                   struct grid ***fields,
                   int nparticles,
                   double dt,
                   int dump);

void resetfield_rho_j(struct grid ***fields, int size);

#endif
