# Known Bugs

Known bugs with `gfortran` compiler:
- (version 8.0), an error is issued ("Fortran runtime error: Bad repeat count in item * of list input"), if
to_file keyword is not last in respective data_array block.
- (version 6.0-8.0), the file mpr_global_paramaters.nml raises an end-of-file-error.
Strangely, it works when adding another line with a backslash at the end of the file.

