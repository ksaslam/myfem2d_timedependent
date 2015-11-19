#include <iostream>
#include <cstdlib>
#include <fstream>
#include <ctime>
#include "constants.hpp"
#include "parameters.hpp"
#include "matprops.hpp"
#include "utils.hpp"
#include "solver.hpp"
#include "fields.hpp"
#define HEIGHT 200
#define WIDTH 200


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


 int n,m,i, iter, maxiter;
 //const int *conn; 
 maxiter=1500;
 double sum;
 
 int number_of_nodes=3;

 double **matrix_first,*vector, *vector_second, *vector_guess;


 vector = (double *)malloc(HEIGHT*sizeof(double));    // this is the x-vector
 vector_second = (double *)malloc(HEIGHT*sizeof(double)); // this is the b-vector
 vector_guess = (double *)malloc(HEIGHT*sizeof(double));	// this is a guess vector that will converge to solution
 
	
 matrix_first=(double **) malloc(HEIGHT*sizeof(double *));  // this is matrix A where Ax = b

 for(i=0;i<HEIGHT;i++)
    matrix_first[i]=(double *) malloc(WIDTH*sizeof(double));  // poplationg matrix A with tridiagnol elements -2,1,1
  
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
    srand(time(NULL));                                 // Populating vector x with random values using random generator
  	vector[n]= ((double)rand() / (double)(RAND_MAX));
  	vector_second[n]= 0. ;
  	//printf("\nThe vector is\n");
  	//printf("%f\t",vector[n]);
  }	
 multiply_matrix (matrix_first, vector, vector_second,WIDTH);   // multipying matrix A and vector x to get vector b.
 
 for (n=0; n<HEIGHT; n++)
  {	
  	vector_guess[n]= vector_second[n] ;                      // giving b vector as a initial guess to test. 
  	
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
 	
  cg_solve(matrix_first, vector_guess,vector_second, WIDTH);  // caling jacobi method to solve for x giving an initial guess.
 }
 
 output << "\nThe vector converged after solver is \n";

 for (n=0;n<HEIGHT;n++)
 {
 	 output << vector_guess[n] << " " ;
 }	 
 for(i=0;i<HEIGHT;i++)
    free(matrix_first[i]);
 free(matrix_first);
 free(vector);
 free(vector_second);
 free(vector_guess);
 
 const int *conn= (*var.connectivity)[1];

  
 output << "\n the number of elements are \n";
 //output << conn[1] << " " ;
 
}





 




