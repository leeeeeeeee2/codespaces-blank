# Configuration

[TOC]

The general concept of the namelist organization is object-oriented. 
There are four different objects: 
[Coordinates](06_namelist.md#coordinates_mpr_nml), 
[Data Arrays](06_namelist.md#data_arrays_mpr_nml), 
[Parameters](06_namelist.md#parameters_mpr_nml) and 
[Upscalers](06_namelist.md#upscalers_mpr_nml).
Each of these objects can be referenced at other locations by name.

Thus, please make sure to not have duplicated names and pay attention to lower/uppercase
so a unique reference can be established.

## Coordinates {#coordinates}

Coordinates are the values associated to an axis of the multi-dimensional data arrays.
Coordinate values must always be always monotonically increasing or decreasing. 
Their name `coord_name` and reference `coord_stagger` must always be defined in the namelist. Apart from that,
they can be initialized in many ways:

1. If `coord_from_file` is specified in the namelist, the values are read from file.
A coordinate variable should not be projected (tilted) and depend on another coordinate variable.
A coordinate can be dependent on 2 dimensions, as long as they are rectangular and
the relevant dimension (of the 2) for this axis is set in the `coordinate_aliases`.
If the coordinate has varying cell widths, then a bound has to be set. This can
either be done in `coord_from_values_bound` or directly in the file using the conventions: 
http://cfconventions.org/Data/cf-conventions/cf-conventions-1.7/cf-conventions.html#cell-boundaries
If no bound is set, the coordinate is not properly initialized.

2. If `coord_from_values` is specified in the namelist, the values are set directly.
If the coordinate has varying cell widths or has only two values (no step size can be inferred), 
then a bound `coord_from_values_bound` has to be set.
If no bound is set, the coordinate is not properly initialized.

3. If `coord_from_range_step` and `coord_from_range_start` and `coord_from_range_count` are specified in the namelist, 
the values are computed directly. Consider `coord_stagger`, when setting this.
If any of these is not set, the coordinate is not properly initialized.

4. If a coordinate is not properly initialized, its missing properties can be inferred
from a source coordinate. A source coordinate is a coordinate which is 
in the same `coordinate_aliases` group, is properly initialized and is a coordinate
of a data array referred to by `from_data_arrays` (this is where the linking happens).
Multiple combinations of missing properties can occur:
   * `coord_from_values_bound` was not provided (e.g. setting from_file): the missing bound is taken from the source coordinate
   * `coord_from_range_start` and `coord_from_range_count` were not provided: the missing bounds are taken from the source coordinate
   and a new count and start value are calculated
   * `coord_from_range_start` and `coord_from_range_step` were not provided: the missing bounds are taken from the source coordinate
   and a new step and start value are calculated

## Transfer functions {#transfer_functions}

#### general form

We support mathematical expressions for transfer function
names in the following form:

`<p|data> <op> <p|data>`

where `<p|data>` denotes either parameter or data. 

#### possible operators

`<op>` denotes one
of the following operators: 
`+` , `-` , `**` , `*` , `/` , `(`, `)`, `exp`, `log10`, `log`, `if`, `else`, `then`, `where`, `end`, `<=` , `<` , `>=` , `>` , `==` , `.and.` , `.or.` , `sin` , `cos` , `tan` , `abs` , `acos`, `max` , `min` , `asin`, `atan`, `atan2`, `cosh`, `sinh`, `tanh`, `sqrt`.
They are directly parsed to Fortran code.

#### parameters

You can provide parameters locally for a given data array. This is recommended for transfer functions that occur once only. An example:

	&Data_Arrays
	  [...]
      name(3) = "porosity"
      transfer_func(3) = "p1 + p2 * sand + p3 * clay"
      from_data_arrays(1:2,3) = "sand", "clay"
      from_parameter_values(1:3,3) = 0.5, 1.2, 0.8
	/
	
You can also reference the global parameters (as defined in the Parameters section of
the namelist configuration files) by a vector of names. This is recommended for transfer functions that occur more than once and have many parameters (and long names). An example:

	&Data_Arrays
	  [...]
      name(3) = "porosity"
      transfer_func(3) = "p1 + p2 * sand + p3 * clay"
      from_data_arrays(1:2,3) = "sand", "clay"
      from_parameter_names(1:3,3) = "PTF1_p1", "PTF1_p2", "PTF1_p3"
	/
	&Parameters
	  parameter_names(1:3) = "PTF1_p1", "PTF1_p2", "PTF1_p3"
	  parameter_values(1:3) = 0.5, 1.2, 0.8
	/

You can also reference the global parameters (as defined in the Parameters section of
the namelist configuration files) directly by name. This is recommended for transfer functions that occur more than once and have not so many parameters. An example:

	&Data_Arrays
	  [...]
      names(3) = "porosity"
      transfer_funcs(3) = "PTF1_p1 + PTF1_p2 * sand + PTF1_p3 * clay"
      from_data_arrays(1:2,3) = "sand", "clay"
	/
	&Parameters
	  parameter_names(1:3) = "PTF1_p1", "PTF1_p2", "PTF1_p3"
	  parameter_values(1:3) = 0.5, 1.2, 0.8
	/

#### conditional clauses

An if-clause needs to be used as following: 

`if ( _logical_expression_ ) then ... else if ( _logical_expression_ ) then ... else ...`

A where-clause needs to be used as following: 

`where ( _logical_expression_ ) ... else where ( _logical_expression_ ) ... else ...`

In both cases, the `else if`/`else where` are optional and also the `else` is optional. 
Please remember to handle all possible cases of an if or where-clause. 
Generally, a where-clause must be used if arrays are used in the logical expression and 
an if-clause if scalars are used in the logical expression.
You need not consider nan-values in the transfer function as the transfer function is applied to nan-free (masked) arrays only.

Fortran requires any where-clause to have an "else where" statement.
If the user omits that, MPR inserts an "else where\nfunc_result(:) = nodata_dp" additionally.

Please make sure that you do not accidentally name your parameters or data like one of the operators
(e.g. good candidates are `max` or `min`). Also do not use parameter names x{}, where {} stands for any number.

Whitespaces can be omitted, but not between words (e.g. `elsewhere` is not valid).
Multi-line transfer function names are also okay.
There must not be any line break character and make sure to begin the next line 
with a whitespace if necessary.

#### missing values

The transfer function is applied to non-masked values only.
You can specify missing values in the input netCDF files as specified by the [CF Conventions](http://cfconventions.org/Data/cf-conventions/cf-conventions-1.8/cf-conventions.html#missing-data).
If a transfer function introduces new missing values (by default set to `-9999.0`) in a data array, this is not accounted for in the default settings. 
In order to include the new missing values in the mask after the application of the transfer function, set the property `check_for_nodatavalue` in the [`main` section](06_namelist.md#main_mpr_nml) to `.true.`.

#### transfer func labels

You can provide a label for a transfer function, it is currently only used in the graphviz visualizaiton tools but not in MPR itself: 

	&Data_Arrays
	  [...]
      name(3) = "porosity"
      transfer_func(3) = "TF1_p1 + TF1_p2 * sand + TF1_p3 * clay"
      transfer_func_label(3) = "MyFineTF_for_porosity"
      from_data_arrays(1:2,3) = "sand", "clay"
	/
	&Parameters
	  parameter_names(1:3) = "TF1_p1", "TF1_p2", "TF1_p3"
	  parameter_values(1:3) = 0.5, 1.2, 0.8
	/

## Special configurations

### Using coordinate values as data arrays

Sometimes it is required to use coordinate values directly as data arrays.
You might for example want to use the depth/height of a cell as input to a transfer function.

	&Main
	  coordinate_aliases(:,1) = "lon", "lon_out", "x"
	  coordinate_aliases(:,2) = "lat", "lat_out", "y"
	  coordinate_aliases(:,3) = "depth", "depth_out", "z"
	/
	&Data_Arrays
	  [...] 
	  name(5) = "data_array5"
	  transfer_func(5) = "data_array1 * (z.upper_bound - z.lower_bound) * 1000.0"
	  from_data_arrays(1:3,5) = "data_array1", "z.lower_bound", "z.upper_bound"
	  target_coord_names(1:3,5) = "lon_out", "lat_out", "depth_out"
	  upscale_ops(1:3,5) = "-1.0", "-1.0", "1.0"
	/

It is possible to use the upper and lower bound of a coordinate.
The syntax to do so is `<coord_name>.upper_bound` and `<coord_name>.lower_bound`.
`coord_name` can refer to any coordinate name specified in the [`main` section](06_namelist.md#main_mpr_nml)
under `coordinate_aliases`. However, the exact values are then taken from the 
input data_array's coordinate associated with `coordinate_aliases`. In this example,
`data_array1` is the only input data array and its `z`-dimension might be named
`depth`. As a result, the coordinate values of `depth` are then used for the transfer function.

They will then be broadcast
to the shape of the other data_array (`data_array1` in the example).
Note that there must always be another data array in `from_data_arrays`, but
it does not need to appear in `transfer_func`.

### Reading data arrays from (unstructured or curvilinear) 2D grids

MPR can read 2D coordinates variables in the SCRIP format ([SCRIP Grid File Format](https://earthsystemmodeling.org/docs/nightly/develop/ESMF_refdoc/node3.html#SECTION03028100000000000000)).
Each coordinate should reside in a different netcdf file and should be initialized like this:

	&Main
	    coordinate_aliases(1:4,1) = 'lon', 'my2Dgrid', 'x', 'west_east'
	    coordinate_aliases(1:4,2) = 'lat', 'my2Dgrid', 'y', 'south_north'
	    (...)
	/

	&Coordinates
	    coord_name(1) = 'my2Dgrid'
	    coord_stagger(1) = 'center'
	    coord_from_file(1) = 'path/my2Dgrid.nc'
	    coord_sub_dims(1:2,1) = 'lat', 'lon'
	/

The following variables are read: `grid_size`, `grid_corners`, `grid_center_lon`', `grid_center_lat`', `grid_corner_lon`', `grid_corner_lat`. 
The units attribute should be common for all variables and is read from the `grid_center_lon` variable or alternatively from the global attribute.
It is important that the netcdf file containing the coordinate needs to have the `title` attribute with the coordinate name.
The `coord_sub_dims` section is important for MPR to associate the correct source dimension when remapping multiple singular (1D) dimensions to 2D grid.


### Remapping polygon based data

#### External weights file
The remapping of polygons is fully supported, if the user provides external files with
the subcell ids and weights (use tools like cdo, SCRIP or ESMF). 
There is an example provided, where the two dimensions
`poly_id` and `hruid` are linked through various arrays containing information on the
remapping. 
It requires a vector on the number of subcells for each cell (`n_subcells_field_name`),
an array of the weights of each subcell for each cell (`weights_field_name`) and
an array of the subcells ids used for each cell (`subcell_ids_field_name`):

	&Main
	  read_weights = .True.
	/
	&Upscalers
	  upscaler_name(1) = "polyid__to__hruid"
	  from_weight_file(1) = "./test/spatialweights_600m_hawaii-1km-land_mod.nc"
	  subcell_ids_field_name(1) = "intersector"
	  weights_field_name(1) = "weight"
	  n_subcells_field_name(1) = "overlaps"
	/

Alternatively, the user can also provide the weights and ids in the [SCRIP-based netcdf file format](https://earthsystemmodeling.org/docs/nightly/develop/ESMF_refdoc/node3.html#SECTION03028100000000000000).
This sparse format needs to have the following variables `src_address`, `dst_address` and  `remap_matrix` (instead of `col`, `row` and `S` in the link).

	&main
	  read_weights = .True.
	/
	&Upscalers
	  upscaler_name(1) = "polyid__to__hruid"
	  from_weight_file(1) = "./test/spatialweights_SCRIP.nc"
	/

#### Simple nearest source to destination algorithm

If the user does not provide a weight file, MPR uses a simple approximation.
It assigns an equal weight to each source cell whose center point lies within the target cell. 
It issues a warning if there is not a sufficiently high number of source cells for a target cell to justify this approximation.
