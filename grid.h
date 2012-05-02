#ifndef GRID_H
#define GRID_H

struct grid; 

typedef struct grid
{
    double E[3];          // electric field,
    double B[3];          // magnetic field,
    double rho;           // charge density and
    double J[3];          // current density.
  }grid;

#endif
