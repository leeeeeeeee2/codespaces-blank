### Getting started

This guide walks you through configuring MPR for a basic example of applying
a transfer function to a set of predictor variables and scale it
to a specified target grid. Each of the variables is a multi-dimensional array with references to multiple dimensions (comparable to NetCDF file format).

At this point, it is important to introduce the naming convention of MPR.
Any multidimensional variable is called _DataArray_.
Each dimension of any DataArray needs to be labelled and have values assigned.
This is done by a _Coordinate_. _Coordinates_ are usually one-dimensional and equal-spaced but do not need to be.
In MPR, _Coordinates_ define intervals and not points for each value in a _DataArray_.
For example, soil information are often provided for soil depths.
A _Coordinate_ would for example explicitly need to define the boundaries of each layer, 
like [[0, 0.5], [0.5, 1.0]], so the remapping or upscaling of values works unambiguously.

MPR requires one configuration file `mpr.nml`. Let's begin with the `mpr.nml` file. 

1. Here we define the name of the output file.
   
        &main
          out_filename = './template_out.nc'
        /

2. Check that all predictor variables share the same projection and have each dimension correctly set up using the [CF conventions](https://cfconventions.org/Data/cf-conventions/cf-conventions-1.10/cf-conventions.html#coordinate-types).

3. Check the coordinate variables of each of the predictor variables.
    Each coordinate variable needs to get an entry 'coordinate_aliases' in the configuration file 
    specifying the dimension they refer to (for example,
    the predictors come in 3D with dimensions `lon`, `lat` and `depth`). 
    Each axis or dimension needs to be defined so MPR knows, which coordinates refer to the same axis during the upscaling process:

        &main
          out_filename = './template_out.nc'
          coordinate_aliases(:,1) = "lon", "x"
          coordinate_aliases(:,2) = "lat", "y"
          coordinate_aliases(:,3) = "depth", "z"
        /
	
4. Now we need to configure all the actual predictor variables.
    They go into an extra block (for example, the variables names
    are `sand` and `clay`).
   
        &Data_Arrays
          name(1) = "sand"
          from_file(1) = "./src/tests/files/sand.nc"
          to_file(1) = .false.
		  
          name(2) = "clay"
          from_file(2) = "./src/tests/files/clay.nc"
          to_file(2) = .false.
        /
   
    By default, all variables are written to the output file.
    If we do not want a predictor variable in the output file,
    we need to add `to_file(...) = .false.` as in the example above.
   
5. Now we can configure our output variable that needs to be calculated by
    a transfer function. Let us assume, we want the transfer function
    `a + b * sand + c * clay` to calculate `porosity`. We add this to the
    `&Data_Arrays` block. 
   
        &Data_Arrays
          [...]
          name(3) = "porosity"
          transfer_func(3) = "a + b * sand + c * clay"
          from_data_arrays(1:2,3) = "sand", "clay"
        /

    We tell MPR to use `sand` and `clay` as predictors by setting them
    at `from_data_arrays`. [Certain conventions](05_configuration.md) apply for formulating a transfer function,
    but generally follow math notations as used in the Fortran language.
   
6. Another thing we now need is to give meaning
    to the parameters `a`, `b` and `c` that were defined in the transfer function during the previous step. 
    They get another block in the configuration file.
   
        &Parameters
          parameter_names(1:3) = "a", "b", "c"
          parameter_values(1:3) = 0.5, 1.2, 0.8
        /
   
7. Now the upscaling step can be defined. We want the `porosity` field 
    be scaled to a certain grid with the coordinates
    `lon_out`, `lat_out` and `depth_out`. For each of the coordinates
    the aggregation shall be based on the arithmetic mean (with 
    power mean parameter `p=1.0`). This information can be added to the Data_Array block.
		
        &Data_Arrays
          [...]
		  
          name(3) = "porosity"
          transfer_func(3) = "a + b * sand + c * clay"
          from_data_arrays(1:2,3) = "sand", "clay"
          target_coord_names(1:3,3) = "depth_out", "lat_out", "lon_out"
          upscale_ops(1:3,3) = "1.0", "1.0", "1.0"
        /

8. MPR now needs to know, how the target coordinates are defined. You can
    configure that in [multiple ways](05_configuration.md). In any case, we have to define
    another block in the configuration file. The coordinates can for
    example be taken from an existing netcdf file containing these 
    coordiante variables. 
  
        &Coordinates
          coord_name(1:3) = "lon_out", "lat_out", "depth_out"
          coord_stagger(1:3) = "center", "center", "end"
          coord_from_file(1:3) = "./dim_example.nc", "./dim_example.nc", "./dim_example.nc"
        /

    Each coordinate name of the target dimensions defined for the `porosity` data array
    must be supplied. `coord_stagger` defines to which reference point in the cell the coordinate
    value is referring (options: `start`, `center` (default) or `end`).
    If there are no reference files for the target grids available, there are multiple other ways available to initialize coordinates.
    Coordinates can be defined by defining a step size (`coord_from_range_step`) (implicitly depends on source coordinate) 
    or concrete values (`coord_from_values` and `coord_from_values_bound`) (explicit definition).  

        &Coordinates
          coord_name(1:3) = "lon_out", "lat_out", "depth_out"
          coord_stagger(1:3) = "center", "center", "end"
          coord_from_range_step(1) = 0.03125
          coord_from_range_step(2) = 0.03125
          coord_from_values_bound(3) = 0.0
          coord_from_values(:,3) = 0.1, 0.3, 0.6, 1.0
        /
    0.03125

9. In order to link the target coordinates to the correct source coordinates, we have to
    alter the first block (`main`).
  
        &main
          out_filename = './template_out.nc'
          coordinate_aliases(:,1) = "lon", "lon_out", "x"
          coordinate_aliases(:,2) = "lat", "lat_out", "y"
          coordinate_aliases(:,3) = "depth", "depth_out", "z"
        /
   
10. Our final namelist file `mpr.nml` now looks like this:

        &main
          out_filename = './template_out.nc'
          coordinate_aliases(:,1) = "lon", "lon_out", "x"
          coordinate_aliases(:,2) = "lat", "lat_out", "y"
          coordinate_aliases(:,3) = "depth", "depth_out", "z"
        /
		
        &Coordinates
          coord_name(1:3) = "lon_out", "lat_out", "depth_out"
          coord_stagger(1:3) = "center", "center", "end"
          coord_from_range_step(1) = 0.03125
          coord_from_range_step(2) = 0.03125
          coord_from_values_bound(3) = 0.0
          coord_from_values(:,3) = 0.1, 0.3, 0.6, 1.0
        /
		
        &Parameters
          parameter_names(1:3) = "a", "b", "c"
          parameter_values(1:3) = 0.5, 1.2, 0.8
        /
		
        &Data_Arrays
          name(1) = "sand"
          from_file(1) = "./src/tests/files/sand.nc"
          to_file(1) = .false.
		  
          name(2) = "clay"
          from_file(2) = "./src/tests/files/clay.nc"
          to_file(2) = .false.
		  
          name(3) = "porosity"
          transfer_func(3) = "a + b * sand + c * clay"
          from_data_arrays(1:2,3) = "sand", "clay"
          target_coord_names(1:3,3) = "depth_out", "lat_out", "lon_out"
          upscale_ops(1:3,3) = "1.0", "1.0", "1.0"
        /
     
    You can now move (parts of) the parameters in the `&Parameters` section to yet another configuration file (e.g. `mpr_global_parameter.nml`),
    e.g. if you would want these to be changed frequently by an optimization algorithm.
    Follow the same format as used in the `&Parameters` section.
    The indices of `parameter_values` and `parameter_names` define the Parameter entry.
    Note that also numeric constants need to be assigned explicitly (e.g. `1.0` for the transfer function `a * (b - 1.0)`).

11. The next step is setting up the environment for MPR. 
    Please refer to the [installation section](03_installation.md).

12. Now that you have configured MPR, running MPR will raise
    an error as the Fortran code does not know how to interpret
    the transfer function as it is a comiled language. 
    In order to add it to the source code, do as
    told by the error message and load a Python environment where
    the library `f90nml` is installed (see [Dependencies section](02_dependencies.md)).
     
    Running `python -m src_python.pre_proc.update_tfs_in_fortran_source` will add the transfer
    functions found in the `mpr.nml` directly to the source code.
    Display the options with `python -m src_python.pre_proc.update_tfs_in_fortran_source --help`.
    Alternatively, pass the path to the namelist directly to cmake by `-DGENERATE_TFS_CONFIG=<path_to_nml>`.
    You can then compile the newly added transfer functions and MPR will be able to execute the
    minimal example.