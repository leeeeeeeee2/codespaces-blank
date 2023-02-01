!> \file mo_mpr_file.f90
!> \brief provides file names and units for MPR
!> \authors robert schweppe
!> \date jan 2018

module mo_mpr_file

  implicit none

  !> default name for main config file (namelist)
  character(len = *), parameter :: filenamenamelistmprdefault = 'mpr.nml'
  !> unit for main config file (namelist)
  integer, parameter :: unamelistmpr = 80

  !> default name for optional parameter config file (namelist)
  character(len = *), parameter :: filenamenamelistparamdefault = 'mpr_global_parameter.nml'
  !> unit for optional parameter config file (namelist)
  integer, parameter :: unamelistparam = 81

end module mo_mpr_file
