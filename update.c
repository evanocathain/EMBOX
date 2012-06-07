#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#include"grid.h"
#include"particles.h"
#include"update.h"

void update_field_current(struct particles *charges,
                          struct grid ***fields,
                          int nparticles,
                          double dx,
                          double dy, 
                          double dz){

    double q = 1.60217646e-19;
    int i, x_pos, y_pos, z_pos;  // loop var

    for(i=0; i<nparticles; i++){
      // get approximate charge position in 3-D
      x_pos=(int)(charges[i].x[0]/dx);
      y_pos=(int)(charges[i].x[1]/dy);
      z_pos=(int)(charges[i].x[2]/dz);

        //fprintf(stdout, "%i\t%i\t%i\n", x_pos, y_pos, z_pos);
      // update rho and J for the relevant field location
      fields[x_pos][y_pos][z_pos].rho += q*charges[i].q;
      fields[x_pos][y_pos][z_pos].J[0]+=q*charges[i].q*charges[i].u[0];
      fields[x_pos][y_pos][z_pos].J[1]+=q*charges[i].q*charges[i].u[1];
      fields[x_pos][y_pos][z_pos].J[2]+=q*charges[i].q*charges[i].u[2];
    }

}

void update_field_strength(struct grid ***fields,
                 int size,
                 double dx, 
                 double dy, 
                 double dz,
                 double dt, 
                 int dump){

    // derivatives of E and B fields
    double ddx_Ey, ddx_Ez;
    double ddy_Ex, ddy_Ez;
    double ddz_Ex, ddz_Ey;
    double ddx_By, ddx_Bz;
    double ddy_Bx, ddy_Bz;
    double ddz_Bx, ddz_By;

    // loop variables
    int i, j, k;

    // permeability of the vacuum in SI units
    const double mu0 = 4.0*M_PI*10.0e-7;


    /* 
     * SPATIAL DERIVATIVES FOR ELECTRIC / MAGNETIC FIELD
     */
    for (i=0; i<size; i++){
        for (j=0; j<size; j++){
            for (k=0; k<size; k++){

                /* X DERIVATIVES */
                if (i==0){
                    ddx_Ey = ddx_Ez = 0.0;
                    ddx_By = ddx_Bz = 0.0;
                }
                else{
                    /* ELECTRIC FIELD*/
                    ddx_Ey = (fields[i][j][k].E[1] - fields[i-1][j][k].E[1])/dx;
                    ddx_Ez = (fields[i][j][k].E[2] - fields[i-1][j][k].E[2])/dx;

                    /* MAGNETIC FIELD*/
                    ddx_By = (fields[i][j][k].B[1] - fields[i-1][j][k].B[1])/dx;
                    ddx_Bz = (fields[i][j][k].B[2] - fields[i-1][j][k].B[2])/dx;
                }

                /* Y DERIVATIVES */
                if (j==0){
                    ddy_Ex = ddy_Ez = 0.0;
                    ddy_Bx = ddy_Bz = 0.0;
                }
                else{
                    /* ELECTRIC FIELD*/
                    ddy_Ex = (fields[i][j][k].E[0] - fields[i][j-1][k].E[0])/dy;
                    ddy_Ez = (fields[i][j][k].E[2] - fields[i][j-1][k].E[2])/dy;

                    /* MAGNETIC FIELD*/
                    ddy_Bx = (fields[i][j][k].B[0] - fields[i][j-1][k].B[0])/dy;
                    ddy_Bz = (fields[i][j][k].B[2] - fields[i][j-1][k].B[2])/dy;
                }

                /* Z DERIVATIVES */
                if (k==0){
                    ddz_Ex = ddz_Ey = 0.0;
                    ddz_Bx = ddz_By = 0.0;
                }
                else{
                    /* ELECTRIC FIELD*/ 
                    ddz_Ex = (fields[i][j][k].E[0] - fields[i][j][k-1].E[0])/dz;
                    ddz_Ey = (fields[i][j][k].E[1] - fields[i][j][k-1].E[1])/dz;

                    /* MAGNETIC FIELD*/
                    ddz_Bx = (fields[i][j][k].B[0] - fields[i][j][k-1].B[0])/dz;
                    ddz_By = (fields[i][j][k].B[1] - fields[i][j][k-1].B[1])/dz;
                }


                /* do the field updates */
                fields[i][j][k].E[0] += dt*(ddy_Bz - ddz_By - mu0*fields[i][j][k].J[0]);
                fields[i][j][k].E[1] += dt*(ddz_Bx - ddx_Bz - mu0*fields[i][j][k].J[1]);
                fields[i][j][k].E[2] += dt*(ddx_By - ddy_Bx - mu0*fields[i][j][k].J[2]);
                fields[i][j][k].B[0] += dt*(ddz_Ey - ddy_Ez);
                fields[i][j][k].B[1] += dt*(ddx_Ez - ddz_Ex);
                fields[i][j][k].B[2] += dt*(ddy_Ex - ddx_Ey);

                if (dump == 1) {
                    printf("%d %d %d %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf\n",i,j,k,fields[i][j][k].E[0],fields[i][j][k].E[1],fields[i][j][k].E[2],fields[i][j][k].B[0],fields[i][j][k].B[1],fields[i][j][k].B[2],fields[i][j][k].J[0],fields[i][j][k].J[1],fields[i][j][k].J[2]);
                } // end of if
            } //end of k loop
        } // end of j loop
    } // end of i loop
}

double trilin_interp_E(struct grid ***fields,
                     struct particles *charges,
                     int size,
                     int unit_vec,
                     double xpos,
                     double ypos,
                     double zpos){

    double c00, c10, c01, c11;
    double c0, c1, c;
    // Using convention from wikipedia
    // http://en.wikipedia.org/wiki/Trilinear_interpolation
    // (don't pretend you would have looked it up anywhere else!)
    double xd, yd, zd;
    int x0, y0, z0;
    int x1, y1, z1;
    double sum;

    x0 = floor(xpos);
    x1 = ceil(xpos);

    y0 = floor(ypos);
    y1 = ceil(ypos);

    z0 = floor(zpos);
    z1 = ceil(zpos);
    //fprintf(stderr, "xp %lf yp %lf zp %lf \n", xpos, ypos, zpos);

    // Need to check if floor == ceil, the code will break.
    // ie is the point on a plane, an edge, or a vertex in the grid?
    // THIS IS HORRIBLE
    if (x0 == x1){
        if (y0 == y1){
            if (z0 == z1){
                // on a vertex, return the value there
                return fields[x0][y0][z0].E[unit_vec];
            }
            else{
                // on an edge - average in z space
                sum = fields[x0][y0][z0].E[unit_vec] + fields[x0][y0][z1].E[unit_vec];
                return sum/2.0;
            }
        }
        else if (z0 == z1){
            // edge - average in y space
            sum = fields[x0][y0][z0].E[unit_vec] + fields[x0][y1][z0].E[unit_vec];
            return sum/2.0;
        }
        else{
            // on a plane - average four points
            sum = fields[x0][y0][z0].E[unit_vec] + fields[x0][y1][z0].E[unit_vec]
                + fields[x0][y0][z1].E[unit_vec] + fields[x0][y1][z1].E[unit_vec];
            return sum/4.0;
        }
    }
    else if (y0 == y1){
        if (z0==z1){
            // return x-space average
            sum = fields[x0][y0][z0].E[unit_vec] + fields[x1][y0][z0].E[unit_vec];
            return sum/2.0;
        }
        else{
            // on a plane - average four points
            sum = fields[x0][y0][z0].E[unit_vec] + fields[x1][y0][z0].E[unit_vec]
                + fields[x0][y0][z1].E[unit_vec] + fields[x1][y0][z1].E[unit_vec];
            return sum/4.0;
        }
    }
    else if (z0 == z1){
        sum = fields[x0][y0][z0].E[unit_vec] + fields[x1][y0][z0].E[unit_vec]
                + fields[x0][y1][z0].E[unit_vec] + fields[x1][y1][z0].E[unit_vec];
        return sum/4.0;
    }
    // Reached here? No matches, do full interpolation
    
    xd = (xpos - x0)/(x1 - x0);
    yd = (ypos - y0)/(y1 - y0);
    zd = (zpos - z0)/(z1 - z0);

    //fprintf(stderr, "x1 %i y1 %i z1 %i x0 %i y0 %i z0 %i\n", x1, y1, z1, x0, y0, z0);
    // start interpolating in x
    c00 = fields[x0][y0][z0].E[unit_vec] * (1.0 - xd) +
                    fields[x1][y0][z0].E[unit_vec] * xd;
    c10 = fields[x0][y1][z0].E[unit_vec] * (1.0 - xd) +
                    fields[x1][y1][z0].E[unit_vec] * xd;
    c01 = fields[x0][y0][z1].E[unit_vec] * (1.0 - xd) +
                    fields[x1][y0][z1].E[unit_vec] * xd;
    c11 = fields[x0][y1][z1].E[unit_vec] * (1.0 - xd) +
                    fields[x1][y1][z1].E[unit_vec] * xd;
    //fprintf(stderr, "x %i y %i z %i xd %lf yd %lf zd %lf\n", x0, y0, z0, xd, yd, zd);

    // interpolate in y
    c0 = c00 * (1.0 - yd) + c10 * yd;
    c1 = c01 * (1.0 - yd) + c11 * yd;

    // interpolate in z
    c = c0 * (1.0 - zd) + c1 * zd;

    //fprintf(stderr, "c00 %lf\tc01 %lf\tc10 %lf\tc11 %lf\n", c00, c01, c10, c11);
    //fprintf(stderr, "c: %lf\t%lf\n", c0, c1);
    return c;
}

void update_charge_posns(struct particles *charges,
			 struct grid ***fields,
			 int nparticles,
			 double dt,
			 double dx,
			 double dy,
			 double dz,
			 int size,
			 int dump,
			 FILE *positions){
  
    int i; // loop vars
    int x_pos, y_pos, z_pos;
    double Ex, Ey, Ez, Bx, By, Bz, Exi; // em fields
    double ax, ay, az; // accelerations
    double x_update, y_update, z_update;
    const double q_to_m=1.75882017e11; 
    
    for (i=0; i<nparticles; i++){
        if ((charges[i].x[0] +1.0*dx >= size*dx) || (charges[i].x[1] +1.*dx >= size*dy) || (charges[i].x[2]+1.*dx >= size*dz) || (charges[i].x[0]<= 0.0) || (charges[i].x[1] <= 0.0) || (charges[i].x[2] <= 0.0)){
            charges[i].x[0] = size*dx*0.5;
            fprintf(stdout, "%lf\t%lf\n", charges[i].x[0] +1.0*dx, size*dx);
            charges[i].x[1] = size*dy*0.5;
            charges[i].x[2] = size*dz*0.5;
            charges[i].u[0] = 0.0;
            charges[i].u[1] = 0.0;
            charges[i].u[2] = 0.0;
      }
    else{
        // Currently using (int) cast
        // probably can do something more complex here
        x_pos = (int)(charges[i].x[0]/dx);
        y_pos = (int)(charges[i].x[1]/dy);
        z_pos = (int)(charges[i].x[2]/dz);
        /*
        Ex=fields[x_pos][y_pos][z_pos].E[0];
        Ey=fields[x_pos][y_pos][z_pos].E[1];
        Ez=fields[x_pos][y_pos][z_pos].E[2];
        */
        Bx=fields[x_pos][y_pos][z_pos].B[0];
        By=fields[x_pos][y_pos][z_pos].B[1];
        Bz=fields[x_pos][y_pos][z_pos].B[2];
        // try using the trilinear interpolation
        Ex = trilin_interp_E(fields, charges, size, 0,
                          charges[i].x[0]/dx, 
                          charges[i].x[1]/dy,
                          charges[i].x[2]/dz);
        //fprintf(stderr, "E: %.9lf\t%.9lf\n", Ex, Exi);
                          
        Ey = trilin_interp_E(fields, charges, size, 1,
                          charges[i].x[0]/dx, 
                          charges[i].x[1]/dy,
                          charges[i].x[2]/dz);
        Ez = trilin_interp_E(fields, charges, size, 2,
                          charges[i].x[0]/dx, 
                          charges[i].x[1]/dy,
                          charges[i].x[2]/dz);
        
        
        /* Calculate accelerations */
        ax = (q_to_m)*charges[i].q*(Ex + charges[i].u[1]*Bz - charges[i].u[2]*By);
        ay = (q_to_m)*charges[i].q*(Ey + charges[i].u[2]*Bx - charges[i].u[0]*Bz);
        az = (q_to_m)*charges[i].q*(Ez + charges[i].u[0]*By - charges[i].u[1]*Bx);
      
        /* update positions & velocities */
        x_update = charges[i].u[0]*dt+(0.5)*ax*dt*dt;
        y_update = charges[i].u[1]*dt+(0.5)*ay*dt*dt;
        z_update = charges[i].u[2]*dt+(0.5)*az*dt*dt;
        //fprintf(stderr,"%lf\t%lf\t%lf\t%lf\n", charges[i].x[0]/dx, charges[i].x[1]/dy, charges[i].x[2]/dz,x_update);
        charges[i].x[0] += x_update;
        charges[i].x[1] += y_update;
        charges[i].x[2] += z_update;
        charges[i].u[0] += ax*dt;
        charges[i].u[1] += ay*dt;
        charges[i].u[2] += az*dt;
        /*      charges[i].u[0] += 0.5*ax*dt;
        charges[i].u[1] += 0.5*ay*dt;
        charges[i].u[2] += 0.5*az*dt;
        charges[i].x[0] += charges[i].u[0]*dt;
        charges[i].x[1] += charges[i].u[1]*dt;
        charges[i].x[2] += charges[i].u[2]*dt;*/

        //      printf("%f %f %f\n",charges[i].x[0],charges[i].x[1],charges[i].x[2]);
    }
    if (dump == 1){
	    fprintf(positions,"%.9lf %.9lf %.9lf\n",charges[i].x[0],charges[i].x[1],charges[i].x[2]);
      }
    }
}

void resetfield_rho_j(struct grid ***fields, int size){
    int i, j, k; // loop vars

    for(i=0; i<size; i++){
      for(j=0; j<size; j++){
	for(k=0; k<size; k++){
	  fields[i][j][k].rho=0.0;
	  fields[i][j][k].J[0]=0.0;
	  fields[i][j][k].J[1]=0.0;
	  fields[i][j][k].J[2]=0.0;
	}
      }
    }
}
