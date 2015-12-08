#include <iostream>
#include <fstream>
#include "constants.hpp"
#include "parameters.hpp"
#include "binaryio.hpp"
#include "fields.hpp"
#include "geometry.hpp"
#include "ic.hpp"
#include "input.hpp"
#include "matprops.hpp"
#include "mesh.hpp"
#include "output.hpp"
#include "rheology.hpp"
#include "solver.hpp"

#ifdef WIN32
#ifdef _MSC_VER
#define snprintf _snprintf
#endif // _MSC_VER
namespace std { using ::snprintf; }
#endif // WIN32

double **matrix_global, *global_forc_vector, *temperature;
void init_var(const Param& param, Variables& var)
{
    var.time = 0;
    var.steps = 0;
}


void init(const Param& param, Variables& var)
{
    std::cout << "Initializing mesh and field data...\n";

    create_new_mesh(param, var);
    create_boundary_flags(var);
    create_boundary_nodes(var);
    create_boundary_facets(var);
    create_support(var);

    allocate_variables(param, var);

    compute_volume(*var.coord, *var.connectivity, *var.volume);
    
    //this function calcluates the shape functions and its derivatives
    compute_shape_fn(*var.coord, *var.connectivity, *var.volume, *var.shpdx,*var.shpdz, *var.shp ,var);

    
    
    // ******************************************* Initializaation some paramters like global stiffness matrix etc ********************
    
    
    global_forc_vector = (double *)malloc(var.nnode*sizeof(double));
    //vector_guess = (double *)malloc(var.nnode*sizeof(double));
    temperature= (double *)malloc(var.nnode*sizeof(double));
    matrix_global=(double **) malloc(var.nnode*sizeof(double *));  // this is matrix A where Ax = b
    for(int i=0;i<var.nnode;i++)
    {
      matrix_global[i]=(double *) malloc(var.nnode*sizeof(double));  // initializing K matrix 
    }  
    std::ofstream output("./globalmatriks.txt");
    //*********************************************************************************************************


    

    initialize_global_matrix (matrix_global, var.nnode);
    
    initialize_global_force_vector(global_forc_vector, var.nnode);
    
    initialize_guess_temperature_vector(temperature, var.nnode); // this is the temparature, it is the temparature array
    

    
      
// ************************************************************************** 

// computing the global stiffness matrix and global force vector. This function returns both values. 
    
    compute_global_matrix( *var.coord, matrix_global, global_forc_vector, *var.shpdx, *var.shpdz, *var.volume, *var.connectivity,*var.shp, var );
   
    get_steady_temperature_cg_solve(matrix_global,temperature, global_forc_vector, var.nnode);
    
    for(int n=0;n<var.nnode;n++)
    {
          (*var.temperature)[n] = temperature[n] ;
           std::cout << "the temparature value is\n";
           std::cout << temperature[n];
           std::cout << "the final temparature value is\n";
           std::cout << (*var.temperature)[n];
    }      


    // temperature should be init'd before stress and strain
    initial_temperature(param, var, *var.temperature);
   // free_global_matrix(matrix_global,var.nnode); 
    free(matrix_global);
    free(global_forc_vector);
}


int main(int argc, const char* argv[])
{
    double start_time = 0;

    //
    // read command line
    //
    if (argc != 2) {
        std::cout << "Usage: " << argv[0] << " config_file\n";
        std::cout << "       " << argv[0] << " -h or --help\n";
        return -1;
    }

    Param param;
    get_input_parameters(argv[1], param);

    //
    // run simulation
    //
    static Variables var; // declared as static to silence valgrind's memory leak detection
    init_var(param, var);

    Output output(param, start_time,
                  (param.sim.is_restarting) ? param.sim.restarting_from_frame : 0);


    init(param, var);


    var.dt = compute_dt(param, var);
    output.write(var, false);

    double starting_time = var.time; // var.time & var.steps might be set in restart()
    double starting_step = var.steps;
    int next_regular_frame = 1;  // excluding frames due to output_during_remeshing

    std::cout << "Starting simulation...\n";
    do {
        var.steps ++;
        var.time += var.dt;

    
        update_temperature(param, var, *var.temperature, *var.ntmp);

        
        if (
            ( (var.steps - starting_step) == next_regular_frame * param.sim.output_step_interval ) ||
            ( (var.time - starting_time) > next_regular_frame * param.sim.output_time_interval_in_yr * YEAR2SEC )
           )
        {
            output.write(var);

            next_regular_frame ++;
        }
    } while (var.steps < param.sim.max_steps && var.time <= param.sim.max_time_in_yr * YEAR2SEC);

    std::cout << "Ending simulation.\n";
    return 0;
}
