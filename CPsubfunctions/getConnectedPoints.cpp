// Function to get the connected points
// Implemented by Sharif, Chowdhury
// Date: 28.11.2008



#include "mex.h"
#include <stdio.h>

void DFS(int i, int j, int  m, int  n, double color, double *B, int d);

void mexFunction(int nlhs,
                 mxArray *plhs[],
                 int nrhs,
                 const mxArray *prhs[])
{
    int m,n,i,j,k;
    double *A, *B, color,k2;
    /* check for arguments here, omitted */
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    A= mxGetPr(prhs[0]);
    /* allocate the answer */
    plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
    B = mxGetPr(plhs[0]);
    k= 0;
   
    for (i= 0; i<n; i++)
    for(j=0;j<m;  j++) 
    {
        *(B +  (i*m+j)) = *(A +  (i*m+j));
    }  
    k2= 0;
    for (i= 0; i<n; i++)
    for(j=0;j<m;  j++) 
    {
        if ( *( B +  (i*m +j) ) < -0.1 )
        {
            color= *(B +  (i*m +j));
          //  printf("k=%lf i=%d j=%d m=%d n=%d data=%lf data2=%lf \n",k2, i,j,m,n, *( B + (i*m +j)), *( B + (i*m +j)+2));
            color = -color;
            *(B +  (i*m+j)) = color;
            DFS(i,j,m,n,color,B, 0);
        }
            k2=k2+1 ;
    }
}


void DFS(int i, int j, int  m, int  n, double color, double *B, int d)
{
    
        double cond= 0;
        //printf("i=%d j=%d m=%d n=%d depth=%d data=%lf \n", i,j,m,n, d,*( B + (i*m+j)));
        if (i<0) return;
        if (j<0) return;
        if (i>=n) return;
        if (j>=m) return;
        cond =  (*( B + (i*m+j))) - color;
        cond= cond*cond;
        if (cond < 0.001 ) 
        {
            d= d+1;
            *( B + (i*m+j))=  0.5; //-1*( *( B + (i*m+j)));
            DFS(i-1,j,m,n,color,B,d);
            DFS(i+1,j,m,n,color,B,d);
            DFS(i,j-1,m,n,color,B,d);
            DFS(i,j+1,m,n,color,B,d);
            DFS(i-1,j+1,m,n,color,B,d);
            DFS(i-1,j-1,m,n,color,B,d);
            DFS(i+1,j-1,m,n,color,B,d);
            DFS(i+1,j+1,m,n,color,B,d);
        }
 }
