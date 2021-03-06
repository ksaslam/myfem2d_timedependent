# Background

This code is a bare-bones version of DynEarthSol3D (https://bitbucket.org/tan2/dynearthsol3d)
to be used in the term project for CERI 7315/8315 Computational Methods for Geodynamics.


# Build

Requirement:
* You will need a recent C++ compiler that supports C++11 standard. (GNU g++
  4.4 or newer version will suffice.)
* You will need a recent version of Boost::Program_options library (1.42 or
  newer version). Instructions for building the library:
  -- Download the source code from www.boost.org
  -- In the untarred source directory, run "./bootstrap.sh"
  -- In the same directory, run "./b2 --with-program_options -q" to build
     the library.
* You will need Python 2.6+ or 3.2+ and the Numpy package.

Build procedure:
* Edit 'Makefile', modify BOOST_ROOT_DIR if you manually built or installed 
  boost library. If you followed the instructions above to build 
  Boost::Program_options library, set BOOST_ROOT_DIR to the untarred boost
  directory.
* Run "make" to build an optimized executable 'myfem2d'.
* Run "make opt=0" to build a debugging executable.


# Run


* Execute "myfem2d inputfile".
* Execute the executable with '-h' flag to see the available input parameters
  and their descriptions.


# Plot

* Run "2vtk.py modelname" to convert the binary output to VTK files.
* Execute 2vtk.py with '-h' flag to see more usage information.
* Some of the simulation outputs might be disabled. Edit 2vtk.py and
  output.cxx to disable/enable them.
* Plot the VTK files with Paraview or LLNL's Visit program.


# Availability

This software, as well as possible updates, is available from the
following URL:
   https://github.com/echoi/myfem2d/


# License

This program is free software: you can redistribute it and/or modify
it under the terms of the MIT / X Windows System license (see the
file LICENSE for the full text).
