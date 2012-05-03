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
