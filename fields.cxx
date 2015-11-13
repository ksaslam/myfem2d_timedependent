#include <iostream>
#include <cstdlib>
#include <fstream>
#include "constants.hpp"
#include "parameters.hpp"
#include "matprops.hpp"
#include "utils.hpp"
#include "solver.hpp"
#include "fields.hpp"


void allocate_variables(const Param &param, Variables& var)
{
    const int n = var.nnode;
    const int e = var.nelem;

    var.volume = new double_vec(e);
    var.temperature = new double_vec(n);

    var.shpdx = new shapefn(e);
    var.shpdz = new shapefn(e);

    var.mat = new MatProps(param, var);
}



void update_temperature(const Param &param, const Variables &var,
                        double_vec &temperature, double_vec &tdot)
{


 const int WIDTH=5;
 const int HEIGHT=5;
 int n,m,i, iter, maxiter;
 maxiter=100;
 double sum;
 double **matrix_first,*vector, *vector_second, *vector_guess;

vector = (double *)malloc(HEIGHT*sizeof(double));
 vector_second = (double *)malloc(HEIGHT*sizeof(double));
 vector_guess = (double *)malloc(HEIGHT*sizeof(double));	
	
 matrix_first=(double **) malloc(HEIGHT*sizeof(double *));

 for(i=0;i<10;i++)
    matrix_first[i]=(double *) malloc(WIDTH*sizeof(double));
  
 for (n=0; n<HEIGHT; n++)
  {	
  	
    for (m=0; m<WIDTH; m++)
    {
       	matrix_first[n][m]=0.;
      if (n==m) {
      		matrix_first[n][m]=-2.;
        }	
	  if (n==m+1) {
      		matrix_first[n][m]=1.;
        }
	  if (n==m-1) {
      		matrix_first[n][m]=1.;
        }		      		      	
    }
  } 

  //printf("\nThe vector is\n");
  for (n=0; n<HEIGHT; n++)
  {	
  	vector[n]= ((double)rand() / (double)(RAND_MAX));
  	vector_second[n]= 0. ;
  	//printf("\nThe vector is\n");
  	//printf("%f\t",vector[n]);
  }	
 multiply_matrix (matrix_first, vector, vector_second,WIDTH);
 
 for (n=0; n<HEIGHT; n++)
  {	
  	vector_guess[n]= vector_second[n] ;
  	
  }	

// Writing the output for testing

 std::ofstream output("./inputmatriks.txt");
 output << "\n The vector guessed is \n";
 for (n=0;n<HEIGHT;n++)
 {
 	 output << vector_guess[n] << " " ;
 }	  

 output << "\nThe Matrix A is \n";
 for (n=0;n<HEIGHT;n++)
 {
	for (m=0;m<WIDTH;m++)
	{
		output << matrix_first[n][m] << " "; 
	}	
	output << "\n";
 }
 output << "\nThe vector x is \n";
 for (n=0;n<HEIGHT;n++)
 {
 	 output << vector[n] << " " ;
 }	 

  for (iter=0; iter< maxiter; iter++)
 {	
 	jordan_solve(matrix_first, vector_guess,vector_second, WIDTH); 
 }
 
 output << "\nThe vector converged after solver is \n";
 for (n=0;n<HEIGHT;n++)
 {
 	 output << vector_guess[n] << " " ;
 }	 

}





 




