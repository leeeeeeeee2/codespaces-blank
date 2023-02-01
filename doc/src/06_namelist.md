## Namelist options

### main section (mpr.nml) {#main_mpr_nml}

- `out_filename` (_String_): output filename. _Required_
- `coordinate_aliases` (_String_): 2-dimensional array with names that refer to the same axis in the coordinate space. This is important for possible upscaling and comparison of coordinates. An example for a row in the array could be: ("x", "lon", "west_east"). _Required_
- `read_weights` (_Boolean_): indicating whether weights for upscaling should be read from file (see [Upscaler](05_configuration.md#upscalers_mpr_nml)). _Optional, default: .false._
- `write_weights` (_Boolean_): indicating whether calculated/read weights for upscaling should be written to a netcdf file in the SCRIP format. _Optional, default: .false._
- `check_for_nodatavalue` (_Boolean_): indicating whether data arrays are checked for newly introduced missing values after the application of the transfer function. _Optional, default: .false._

### Data_Arrays section (mpr.nml) {#data_arrays_mpr_nml}

All following entries are vectors and have to be specified with the index of the data array they belong to (e.g. ```name(2) = 'foo'``` and ```from_data_arrays(1:2,2) = 'bar', 'baz'``` refers to data array 2).
The indices do not have to be consecutive, duplications are not allowed.  

- `name` (_String_): Either the variable name in the netcdf file (if `from_file` is specified), otherwise the target name that will be created for the data array. _Required_
- `from_file` (_String_): netcdf file name that an input array should be read from. _Optional (Either `from_file` or `from_data_arrays` need to be set)_
- `to_file` (_Boolean_): whether data array should be written to `out_filename`. _Optional, default: .true._
- `transfer_func` (_String_): Either name of transfer function to be used or mathematical expression (see conventions below). The transfer functions are contained in the source file ./src/mo_mpr_transfer_func.f90._Optional, default: "identity"_
- `transfer_func_label` (_String_): Label for the transfer function to be used. This property is currently not used other than labelling graphviz plots. _Optional_
- `from_data_arrays` (vector of _String_): Vector of data_array names (see `name` entry in [Data arrays](05_configuration.md#data_arrays_mpr_nml)), on which to apply the function (`transfer_func`). _Optional (Either `from_file` or `from_data_arrays` need to be set)_
- `from_parameter_names` (vector of _String_): Vector of parameter names (see `parameter_names` entry in [Parameters](05_configuration.md#parameters_mpr_nml)), to be used with the function (`transfer_func`). _Must be provided, if transfer function uses "global" parameter names (see [Transfer functions](05_configuration.md#transfer_functions)). Optional_
- `from_parameter_values` (vector of _Float_): Vector of parameter values to be used with the function (`transfer_func`). _Must be provided, if transfer function uses "local" parameter values (see [Transfer functions](05_configuration.md#transfer_functions)). Optional_
- `target_coord_names` (vector of _String_): Vector of target coordinates (see `coord_name` entry in [Coordinates](05_configuration.md#coordinates_mpr_nml)), if upscaling is desired. _Optional_
- `upscale_ops` (vector of _String_): Vector of names of upscaling operator (e.g., `"sum"`, `"var"`, `"std"`, `"max"`, `"min"`) or value of p-norm specified as string (e.g., `"-1"`, `"0"`, `"1.574"`) for each item in `target_coord_names`. _Required, if `target_coord_names` is set. Optional_
- `limits` (_Float_): Vector of limits (min, max) to apply to values after (optional) transfer function application. Both, only one or no limits can be supplied, defaults to no limit. _Optional_

### Coordinates section (mpr.nml) {#coordinates_mpr_nml}

All following entries are vectors and have to be specified with the index of the data array they belong to (e.g. `coord_name(2) = 'foo'` and `coord_from_values(:,2) = 0.0,1.5,4.0`refers to Coordinate 2).
The indices do not have to be consecutive, duplications are not allowed.
An overview on how to define a coordinate can be found [here](05_configuration.md#coordinates). 

- `coord_name` (_String_): Name of the coordinate. _Required_
- `coord_stagger` (_String_): Point of cell the coordinate values refer to, either `"start"`, `"center"`, `"end"`. _Optional, default: `"center"`_
- `coord_from_file` (_String_): Filenames of netcdf files that contain the coordinate variable. _Optional_
- `coord_from_values` (_Float_): 1-dimensional array of target coordinate values. Typically, this array should be specified if the target coordinate is irregular. _Optional_
- `coord_from_range_step` (_Float_): Step size for coordinate values. Typically, this value should be specified if the target coordinate is regularly spaced. For a proper initialization, `coord_from_range_start` and `coord_from_range_count` also need to be supplied. _Optional_
- `coord_from_range_start` (_Float_): Value for first cell of coordinate, consider `coord_stagger` when setting. Only needed when `coord_from_range_count` and `coord_from_range_step` are also given (initialization by range). _Optional_
- `coord_from_range_count` (_Integer_): Number of `coord_from_range_step` (not required if `coord_from_values` is given). _Optional_
- `coord_from_values_bound` (_Float_): Outer bound of coordinate values. These are interpreted relative to `coord_stagger`. If `coord_stagger` is `"start"`, then `coord_from_values_bound` denotes the upper bound of the last cell. If `coord_stagger` is `"end"` or `"center"`, then `coord_from_values_bound` denotes the lower bound of the first cell. _Optional_
- `coord_attribute_names` (_String_): 1-dimensional array of coordinate attribute names. They will be written to the coordinate attributes of the NetCDF output file (e.g. `"standard_name"`, `"long_name"`, `"units"`, `"axis"`). The length needs to be equal to `coord_attribute_values`. _Optional_
- `coord_attribute_values` (_String_): 1-dimensional array of coordinate attribute values.  They will be written to the coordinate attributes of the NetCDF output file (e.g. `"depth"`, `"positive downwards upper boundary of soil horizons"`, `"m"`, `"Z"`). The length needs to be equal to `coord_attribute_names`. _Optional_
- `coord_proj_string` (_String_): The string will define the coordinate projection of that coordinate. Currently this property is not used in MPR. (e.g. ```'+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs84=0,0,0'```). _Optional_
- `coord_unit` (_String_): The string will define the unit of the coordinate. (e.g. `"degrees"`, `"radians"` or `"m"`). _Optional_
- `coord_sub_dims` (_String_): 1-dimensional array of sub dimensions for a coordinate. (e.g. `"lat"`, `"lon"`). _Optional_

### Upscalers section (mpr.nml) {#upscalers_mpr_nml}

This section is only required if `read_weights` is `.true.`. 
All following entries are vectors and have to be specified with the index of the data array they belong to (e.g. `upscaler_name(2) = 'foo'` refers to upscaler 2).
The indices do not have to be consecutive, duplications are not allowed.
An example of an upscaler weights file can be found in 
`./src/tests/files/spatialweights_600m_hawaii-1km-land_mod.nc`.

- `upscaler_name` (_String_): specify which source and target coordinate this upscaler is specified for. For example, if source coordinate name is `"soil_pedons"` and target coordinate name is `"hru_id"`, then `"soil_pedons__to__hru_id"`. Note the two underscores around the *to*.
- `from_weight_file` (_String_): netcdf file name that contains the weights.
- `subcell_ids_field_name` (_String_): name of variable in `from_weight_file` that contains the cell ids of the source coordinate for each cell in the target coordinate.
- `weights_field_name` (_String_): name of the variable containing the weights for each cell in the source coordinate for each cell in the target coordinate.
- `n_subcells_field_name` (_String_): name of variable specifying the number of cells in the source coordinate for each cell in the target coordinate.

### Parameters (mpr.nml and/or mpr_global_parameters.nml) {#parameters_mpr_nml}

This section can either be placed in the file _mpr.nml_
or in _mpr_global_parameters.nml_ or both. The idea is here to provide "constant" parameters in the former and
"variable" parameters (e.g. generated by an auto-calibration software) in the latter.

- `parameter_names` (_String_): Array of parameter names.
- `parameter_values` (_Float_): Array of respective parameter values.


