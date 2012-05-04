#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#include"grid.h"

double pick_rand(void)
{
  double num;
  num=(double)rand();
  num=num/(RAND_MAX+0.0);
  return(num); // a number between 0 and 1
}

void free_grid(struct grid ***data, size_t xlen, size_t ylen)
{
    size_t i, j;
    for (i=0; i<xlen; i++){
        if (data[i] != NULL){
            for (j=0; j<ylen; j++){
                free(data[i][j]);
            }
            free(data[i]);
        }
    }
    free(data);
}
