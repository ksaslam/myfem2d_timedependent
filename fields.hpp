#ifndef DYNEARTHSOL3D_FIELDS_HPP
#define DYNEARTHSOL3D_FIELDS_HPP


void allocate_variables(const Param &param, Variables& var);
void reallocate_variables(const Param &param, Variables& var);
void update_temperature(const Param &param, const Variables &var,
                        double_vec &temperature, double_vec &tdot);
void update_strain_rate(const Variables& var, tensor_t& strain_rate);
void update_force(const Param& param, const Variables& var, array_t& force);
void update_velocity(const Variables& var, array_t& vel);
void update_coordinate(const Variables& var, array_t& coord);
void rotate_stress(const Variables &var, tensor_t &stress, tensor_t &strain);
void initialize_global_matrix (double ** matrix_global, int num_nodes);
void free_global_matrix (double ** matrix_global, int num_nodes);
void initialize_global_force_vector(double * global_forc_vector, int num_nodes);
void initialize_guess_temperature_vector(double * guess_vector, int num_nodes);

// added functions 
//double multiply_matrix(double **matrix_first, double *vector_first, double *vector_second, int nelem_in_array); 


#endif
