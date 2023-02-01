!> \file mo_mpr_global_variables.f90

!> \brief Global variables for mpr only

!> \details

!> \authors Robert Schweppe
!> \date Dec 2017

module mo_mpr_global_variables

  implicit none

  private
  ! ------------------------------------------------------------------
  ! Global Containers
  ! ------------------------------------------------------------------
  logical, public :: CHECK_FOR_NODATAVALUE
  logical, public :: WRITE_WEIGHTS
  character(256), public :: OUT_FILENAME

end module mo_mpr_global_variables
