#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>

/* Blest Botttom for discrete data*/
  void estdep(double *x, double *y, int *n, double *blest, double *std)

   {
       int i,j;
       double x0,y0,sum2, sum3, a1, a2, b1, b2, c1, c2, s1, s2, sum11, sum12;
      

       sum3 = 0.0;
       s1 = 0.0;
       s2 = 0.0;
       for(i=0;i<n[0];i++)
       {
           sum11 = 0.0;
           sum12 = 0.0;
           sum2=0.0;
           x0 = x[i]; y0=y[i];

           for(j=0;j<n[0];j++)
           {
               a1 = ( x[j] <= x0);
               a2= ( x[j] < x0);
               b1 = ( y[j] <= y0);
               b2= ( y[j] < y0);
               sum11 += a1/ (double) n[0];
               sum12 += a2/ (double) n[0];
               sum2+=(b1+b2)/ (double) n[0];

           }
          c1 = (sum11)*(sum11)+(sum12)*(sum12)+(sum11)*(sum12) -1.0;
          c2 = sum2 -1.0;
           sum3 += c1*c2;
           s1 +=  c1*c1;
           s2 +=  c2*c2;
       }
        s1 = s1/((double)n[0]);
        s2 = s2/((double)n[0]);
     std[0] = sqrt(s1*s2);
     blest[0] = sum3/((double)n[0]) /std[0];
     }



