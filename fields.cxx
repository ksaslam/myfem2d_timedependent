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
#define HEIGHT 10
#define WIDTH 10


void allocate_variables(const Param &param, Variables& var)
{
    const int n = var.nnode;
    const int e = var.nelem;

    var.volume = new double_vec(e);
    var.temperature = new double_vec(n);
    var.shpdx = new shapefn(e);
    var.shpdz = new shapefn(e);
    var.shp = new shapefn(e);
    //var.matrix_global= new shapefn(n);
    //var.global_forc_vector = new double_vec(n);

    //var.shp2dz = new shapefn(e);
    //var.shp3dx = new shapefn(e);
    //var.shp3dz = new shapefn(e);

    //var.shpdx = new shapefn(e);
    //var.shpdz = new shapefn(e);
    //var.shpdx = new shapefn(e);
    //var.shpz = new shapefn(e);
    //var.shpx = new shapefn(e);
    var.mat = new MatProps(param, var);

}



void update_temperature(const Param &param, const Variables &var,
                        double_vec &temperature, double_vec &tdot)
{


 int n,m,i, iter, maxiter;
 //const int *conn; 
 maxiter=500;
 double sum;
 
 int number_of_nodes=3;

 double **matrix_first, *vector, *vector_second, *vector_guess;
 

 vector = (double *)malloc(HEIGHT*sizeof(double));    // this is the x-vector
 vector_second = (double *)malloc(HEIGHT*sizeof(double)); // this is the b-vector
 vector_guess = (double *)malloc(HEIGHT*sizeof(double));	// this is a guess vector that will converge to solution
 
// ******************************************* test case ********************
 //globforc_vector = (double *)malloc(var.nnode*sizeof(double));
 //global_force_vector= (double *)malloc(var.nnode*sizeof(double));
 //matrix_global=(double **) malloc(var.nnode*sizeof(double *));  // this is matrix A where Ax = b
  //for(i=0;i<var.nnode;i++)
  //{
    //  matrix_global[i]=(double *) malloc(var.nnode*sizeof(double));  // initializing K matrix 
  //}  
// **************************************************************************  
	
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
    //srand(time(NULL));                                 // Populating vector x with random values using random generator
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

  //for (iter=0; iter< maxiter; iter++)
// {	
 	
  cg_solve(matrix_first, vector_guess,vector_second, WIDTH);  // caling jacobi method to solve for x giving an initial guess.
// }
 
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
 
 
 //const int *conn= (*var.connectivity)[1];

  
 output << "\n the number of elements are \n";
 output << *var.connectivity[0][0] << " " ;
 
}

void initialize_global_matrix (double ** matrix_global, int num_nodes)
{
  int i,j;
  

  // matrix_global =(double **) malloc(num_nodes*sizeof(double *));  // this is matrix A where Ax = b


  // for(int i=0;i<num_nodes;i++)
  //   {

  //     matrix_global[i]=(double *) malloc(num_nodes*sizeof(double));  // initializing K matrix 

  //   }  

  for(i=0;i<num_nodes;i++) // initializing the global k matrix
    {       
       for(j=0;j<num_nodes;j++)
       {
         matrix_global[i][j]= 0.; 
       }

    }

}

void initialize_global_force_vector(double * global_forc_vector, int num_nodes)
{
  int j;      
  for(j=0;j<num_nodes;j++)
    {
        global_forc_vector[j]= 0.; 
    }

  
}

void initialize_guess_temperature_vector(double * guess_vector, int num_nodes)
{
  int j;      
  for(j=0;j<num_nodes;j++)
    {
        guess_vector[j]= 0.5; 
    }

  
}


void free_global_matrix (double ** matrix_global, int num_nodes)

{
  int i;
  
  for(i=0;i<num_nodes;i++)
    {
      free(matrix_global[i]);
    }  

  free(matrix_global);
}
 




