#ifndef UPDATE_H
#define UPDATE_H

void update_field_current(struct particles *charges,
                          struct grid ***fields,
                          int nparticles,
                          double dx,
                          double dy,
                          double dz);

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
			 double dx,
			 double dy,
			 double dz,
			 int size,
			 int dump,
			 FILE *positions);

void resetfield_rho_j(struct grid ***fields, int size);

double trilin_interp_E(struct grid ***fields,
                     struct particles *charges,
                     int size,
                     int unit_vec,
                     double xpos,
                     double ypos,
                     double zpos);
#endif
