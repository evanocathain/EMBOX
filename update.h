#ifndef UPDATE_H
#define UPDATE_H

void update_field_current(struct particles *charges,
                          struct grid ***fields,
                          int nparticles);

void update_field_strength(struct grid ***fields,
                int size, 
                double dx, 
                double dy, 
                double dz, 
                double dt,
                int dump);

void update_charge_posns(struct particles *charges, 
                struct grid ***fields,
                int nparticles,
                double dt,
		int dump,
		FILE *positions);

void resetfield_rho_j(struct grid ***fields, int size);

#endif
