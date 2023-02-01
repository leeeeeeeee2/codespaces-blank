! ---------------------------------------------------------------------------------------
!
! MPR program
!
! author: Robert Schweppe
!
! created: January 2018
!
! ----------------------------------------------------------------------------------------
#include "flogging.h"
program mpr

  USE mo_kind, ONLY: i4

  use mo_mpr_data_array, only : MPR_DATA_ARRAYS, delete_parents
  use mo_mpr_global_variables, only : OUT_FILENAME
  use mo_mpr_file, only : filenameNamelistMprDefault, uNamelistMpr, filenameNamelistParamDefault, uNamelistParam
  use mo_mpr_read_config, only : mpr_read_config
  use mo_netcdf, only : NcDataset
  use mo_mpr_constants, only : maxNoDataArrays, maxNameLength
  use mo_mpr_reset, only: reset
  use mo_mpr_reorder_data_array, only: reorder_data_arrays
  use mo_mpr_command_line_args, only: check_mpr_command_line_args
  use flogging

  implicit none

  integer(i4) :: i, j
  character(maxNameLength), allocatable :: deleteFieldNames_local(:)
  character(maxNameLength) :: filenameNamelistMpr, filenameNamelistParam
  type(NcDataset) :: nc

  !call log_disable_cli_arguments
  call log_set_output_time(.true.)
  ! fatal:1, error:2, warn:3, info:4, debug:5, trace:6, subtrace:7
  ! minimum_loglevel = 4

  ! set default values to filepaths and parse command line arguments
  filenameNamelistMpr = filenameNamelistMprDefault
  filenameNamelistParam = filenameNamelistParamDefault
  call check_mpr_command_line_args(filenameNamelistMpr, filenameNamelistParam)

  call mpr_read_config(filenameNamelistMpr, uNamelistMpr, filenameNamelistParam, uNamelistParam)

  ! reorder data array to minimize storage requirement and assure that all required input fields have been read before
  ! call reorder_data_arrays(MPR_DATA_ARRAYS)

  ! open the output file
  nc = NcDataset(OUT_FILENAME, "w")

  do i = 1, maxNoDataArrays
    ! stop if no name is initialized
    if (.not. MPR_DATA_ARRAYS(i)%is_initialized) then
      cycle
    end if

    log_info(*) 'Working on data array: ', trim(MPR_DATA_ARRAYS(i)%name)
    ! execute MPR on each DataArray
    call MPR_DATA_ARRAYS(i)%execute()
    ! write the current parameter to the output file
    if (MPR_DATA_ARRAYS(i)%toFile) then
      call MPR_DATA_ARRAYS(i)%write(nc)
    end if

    ! clean memory
    call delete_parents(MPR_DATA_ARRAYS(i))

  end do

  ! close the output file
  call nc%close()

  ! clean up
  call reset()
  ! ----------------------------
 
  log_info(*) repeat('-',70)
  log_info(*) 'MPR: Finished!'

end program mpr


