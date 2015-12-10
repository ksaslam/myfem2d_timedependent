#ifndef DYNEARTHSOL3D_GEOMETRY_HPP
#define DYNEARTHSOL3D_GEOMETRY_HPP

double dist2(const double* a, const double* b);
void compute_volume(const array_t &coord, const conn_t &connectivity,
                    double_vec &volume);

double compute_dt(const Param& param, const Variables& var);

void compute_shape_fn(const array_t &coord, const conn_t &connectivity,
                      const double_vec &volume,
                      shapefn &shpdx, shapefn &shpdz, shapefn &shp, const Variables& var);

void compute_global_matrix(const array_t &coord, double **matrix_global,double **matrix_mass, double* global_forc_vector,shapefn &shpdx, shapefn &shpdz, const double_vec &volume, const conn_t &connectivity, shapefn &shp, const Variables& var );
//void compute_global_matrix( const array_t &coord, double** matrix_global, double* global_forc_vector,shapefn &shpdx, shapefn &shpdz, const double_vec &volume, const conn_t &connectivity, shapefn &shp, const Variables& var );
void initialize_local_matrix (double ** matrix_local, int num_nodes);
void initialize_local_force_vector(double * local_forc_vector, int num_nodes);
double boundary_value( int boundary_flag);
void forcing_source( double* forc_vector, double * force_node_coordinate,int number_of_nodes);
void initial_temperature_values( double* temperature,int number_of_nodes);
void temperature_at_new_time( double *temperature_new,double **matrix_global, double *global_forc_vector,double **global_mass,double time_step,double *initial_temperature, int number_of_nodes);         



#endif
