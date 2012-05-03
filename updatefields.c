#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include"grid.h"
#include"updatefields.h"

void updatefield(struct grid ***fields,
                 int size,
                 double dx, 
                 double dy, 
                 double dz,
                 double dt, 
                 double mu0,
                 int dump){

    double ddx_Ey, ddx_Ez;
    double ddy_Ex, ddy_Ez;
    double ddz_Ex, ddz_Ey;
    double ddx_By, ddx_Bz;
    double ddy_Bx, ddy_Bz;
    double ddz_Bx, ddz_By;

    int i, j, k;


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
