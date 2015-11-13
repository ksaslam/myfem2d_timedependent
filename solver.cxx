#include <iostream>
#include <cstdlib>
#include "solver.hpp"

double jordan_solve(double **matrix_first, double *vector_guess,double* vector_second, int nelem_in_array) 
//void multiply_matrix(double **matrix_first, double *matrix_second, int nelem_in_array) 
{
  
 int indx,k;
 
 for (indx = 0; indx < nelem_in_array; indx++)
 {
    double product_sum= 0.;
    for (k = 0; k < nelem_in_array; k++) 
    {
        if(k!=indx)
        {    
            product_sum = product_sum + matrix_first[indx][k]*vector_guess[k];
        }    
    }
    vector_guess[indx] =  1./ matrix_first[indx][indx] * ( vector_second[indx]- product_sum); 
 }

return *vector_guess; 
} 


double gauss_solve(double **matrix_first, double *vector_guess,double* vector_second, int nelem_in_array) 
//void multiply_matrix(double **matrix_first, double *matrix_second, int nelem_in_array) 
{
  
 int indx,k;
 
 for (indx = 0; indx < nelem_in_array; indx++)
 {
    double product_sum= 0.;
    double product_sum2=0.;
    for (k = 0; k < nelem_in_array; k++) 
    {
        if(k<=indx-1)
        {    
            product_sum = product_sum + matrix_first[indx][k]*vector_guess[k];
        }    
        if(k>=indx+1)
        {    
            product_sum2 = product_sum2 + matrix_first[indx][k]*vector_guess[k];
        }           
    }
    vector_guess[indx] =  1./ matrix_first[indx][indx] * ( vector_second[indx]- product_sum- product_sum2); 
 }

return *vector_guess; 
} 


double cg_solve(double **matrix_first, double *vector_guess,double* vector_second, int nelem_in_array) 
//void multiply_matrix(double **matrix_first, double *matrix_second, int nelem_in_array) 
{
 double *vector_r, *vector_p, *vector_cof,alpha_p, *vector_bk, b_k_coeffic;
 int indx,k;
 vector_r = (double *)malloc(nelem_in_array*sizeof(double));
 vector_p = (double *)malloc(nelem_in_array*sizeof(double));
 vector_cof = (double *)malloc(nelem_in_array*sizeof(double));
 vector_bk = (double *)malloc(nelem_in_array*sizeof(double));
 multiply_matrix (matrix_first, vector_guess,vector_r, nelem_in_array);
 
 for (k = 0; k < nelem_in_array; k++)
 {
    vector_r[k]= vector_second[k] - vector_r[k];
    vector_p[k]= vector_r[k];
 }  
 
 multiply_matrix(matrix_first, vector_r, vector_cof,nelem_in_array ) ;
 alpha_p = multiply_vector_transpose(vector_r,vector_r,nelem_in_array) / multiply_vector_transpose(vector_r,vector_cof, nelem_in_array) ;

 for (k = 0; k < nelem_in_array; k++)
 {
    vector_guess[k] = vector_guess[k] + alpha_p * vector_p[k];
    multiply_matrix(matrix_first, vector_p, vector_cof,nelem_in_array );
    vector_r[k] = vector_r[k] - alpha_p * vector_cof[k]   ;
    multiply_matrix(matrix_first, vector_p, vector_bk,nelem_in_array ) ;
    b_k_coeffic = multiply_vector_transpose(vector_p,vector_r,nelem_in_array) / multiply_vector_transpose(vector_p,vector_bk, nelem_in_array) ;
 } 
 free(vector_r);
 free(vector_p);
 free(vector_cof);
 free(vector_bk);

}

double multiply_vector_transpose(double *vector_first, double *vector_second, int nelem_in_array) 
//void multiply_matrix(double **matrix_first, double *matrix_second, int nelem_in_array) 
{
  
 int k;
 double result_vector;
 double sum = 0.;
    for (k = 0; k < nelem_in_array; k++) 
    {
        sum = sum + vector_first[k]*vector_second[k];
    }
    result_vector = sum;
 
return result_vector;
}



double multiply_matrix(double **matrix_first, double *vector_first, double *vector_second, int nelem_in_array) 
//void multiply_matrix(double **matrix_first, double *matrix_second, int nelem_in_array) 
{
  
 int indx,k;
 
 for (indx = 0; indx < nelem_in_array; indx++)
 {
    double sum = 0;
    for (k = 0; k < nelem_in_array; k++) 
    {
        sum = sum + matrix_first[indx][k]*vector_first[k];
    }
    vector_second[indx] = sum;
 }
return *vector_second;
} 



