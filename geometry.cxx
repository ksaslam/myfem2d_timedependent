#include <cmath>
#include <limits>
#include <iostream>
#ifdef USE_OMP
#include "omp.h"
#endif

#include "constants.hpp"
#include "parameters.hpp"
#include "matprops.hpp"
#include "utils.hpp"
#include "geometry.hpp"


/* Given two points, returns the distance^2 */
double dist2(const double* a, const double* b)
{
    double sum = 0;;
    for (int i=0; i<NDIMS; ++i) {
        double d = b[i] - a[i];
        sum += d * d;
    }
    return sum;
}


/* Given two points, returns the area of the enclosed triangle */
static double triangle_area(const double *a,
                            const double *b,
                            const double *c)
{
    double ab0, ab1, ac0, ac1;

    // ab: vector from a to b
    ab0 = b[0] - a[0];
    ab1 = b[1] - a[1];
    // ac: vector from a to c
    ac0 = c[0] - a[0];
    ac1 = c[1] - a[1];

#ifndef THREED
    // area = norm(cross product of ab and ac) / 2
    return std::fabs(ab0*ac1 - ab1*ac0) / 2;
#else
    double ab2, ac2;
    ab2 = b[2] - a[2];
    ac2 = c[2] - a[2];

    // vector components of ab x ac
    double d0, d1, d2;
    d0 = ab1*ac2 - ab2*ac1;
    d1 = ab2*ac0 - ab0*ac2;
    d2 = ab0*ac1 - ab1*ac0;

    // area = norm(cross product of ab and ac) / 2
    return std::sqrt(d0*d0 + d1*d1 + d2*d2) / 2;
#endif
}


void compute_volume(const array_t &coord, const conn_t &connectivity,
                    double_vec &volume)
{
    for (std::size_t e=0; e<volume.size(); ++e) {
        int n0 = connectivity[e][0];
        int n1 = connectivity[e][1];
        int n2 = connectivity[e][2];

        const double *a = coord[n0];
        const double *b = coord[n1];
        const double *c = coord[n2];

        volume[e] = triangle_area(a, b, c);
    }
}


double compute_dt(const Param& param, const Variables& var)
{
    // constant dt
    if (param.control.fixed_dt != 0) return param.control.fixed_dt;

    // dynamic dt
    const int nelem = var.nelem;
    const conn_t& connectivity = *var.connectivity;
    const array_t& coord = *var.coord;
    const double_vec& volume = *var.volume;

    double dt_diffusion = std::numeric_limits<double>::max();
    double minl = std::numeric_limits<double>::max();

    for (int e=0; e<nelem; ++e) {
        int n0 = connectivity[e][0];
        int n1 = connectivity[e][1];
        int n2 = connectivity[e][2];

        const double *a = coord[n0];
        const double *b = coord[n1];
        const double *c = coord[n2];

        // min height of this element
        double minh;
        {
            // max edge length of this triangle
            double maxl = std::sqrt(std::max(std::max(dist2(a, b),
                                                      dist2(b, c)),
                                             dist2(a, c)));
            minh = 2 * volume[e] / maxl;
        }

        dt_diffusion = std::min(dt_diffusion,
                                0.5 * minh * minh / var.mat->therm_diff_max);
    }
    
    double dt = dt_diffusion * param.control.dt_fraction;

    if (dt <= 0) {
        std::cerr << "Error: dt <= 0!  " << dt_diffusion << "\n";
        std::exit(11);
    }
    return dt;
}


void compute_shape_fn(const array_t &coord, const conn_t &connectivity,
                      const double_vec &volume, 
                      shapefn &shpdx, shapefn &shpdz, shapefn &shp, const Variables& var)
{
   int e ;
   int node0,node1,node2;
   //const double *coordinates_node0, *coordinates_node1, *coordinates_node2; // To be completed.
   const double beri_centre_x = 1/3;
   const double beri_cente_y = 1/3;
   const double weight= 0.5;
   

   for (e=0; e< var.nelem; e++)
   { 
      shp[e][0]= weight * volume[e]* (1.- beri_centre_x- beri_cente_y ) ;  // this is my phi 1 function
      shp[e][1]= weight* volume[e]* (beri_centre_x);       // phi 2 function
      shp[e][2]= weight* volume[e]* (beri_cente_y);      // phi 3 function

      shpdx[e][0] = -1.;
      shpdz[e][0] = -1.;
      shpdx[e][1] = 1.;
      shpdz[e][1] = 0.;
      shpdx[e][2] = 0.;
      shpdz[e][2] = 1.;
   }

} 
  void compute_global_matrix( double **matrix_global, double *global_forc_vector,shapefn &shpdx, shapefn &shpdz, const double_vec &volume, const conn_t &connectivity, shapefn &shp, const Variables& var )

  {
    int e,i,j;
    double k, b, *forc_vector;
    const int number_of_nodes=3;
    const double weight= 0.5;
    //double **matrix_first;
    forc_vector = (double *)malloc(var.nnode*sizeof(double));


 
    for(i=0;i<var.nnode;i++) // initializing the global k matrix
    {  
       forc_vector[i]= 0.;   
       for(j=0;j<var.nnode;j++)
       {
         matrix_global[i][j]= 0.; 
       }

    }
    k=0.;   // initializing k 
      
    for(e=0;e<var.nelem;i++)
    { 
      for(i=0;i<number_of_nodes;i++)
      {

      }
    
    }

    for(e=0;e<var.nelem;i++)
    {
        for(i=0;i<number_of_nodes;i++)
        {
            b =weight* (shp[e][i] ) * forc_vector[connectivity[e][i]] ; 
            
            for(j=0;j<number_of_nodes;j++)
             {
               k = weight * volume[e]* ( shpdx[e][i] * shpdx[e][j] + shpdz[e][i] * shpdz[e][j] );

               matrix_global[connectivity[e][i]][connectivity[e][j]] += k; 
             }

            global_forc_vector[connectivity[e][i]]+= b;
        }
        
    }
   

             

          


  free(forc_vector);
  } //std::cout << "These ARE  THE coordinates value, test...\n";

   ///std::cout <<  b[1] ;


