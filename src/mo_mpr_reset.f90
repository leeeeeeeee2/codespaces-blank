!> \file mo_mpr_reset.f90

!> \brief reset all global lists of MPR types

!> \details reset all global lists of MPR types
!
!> \authors Robert Schweppe
!> \date Jan 2018

module mo_mpr_reset

  use mo_mpr_util_types, only: MPR_COORD_ALIAS
  use mo_mpr_data_array, only: MPR_DATA_ARRAYS
  use mo_mpr_coordinate, only: MPR_COORDINATES
  use mo_mpr_data_array_upscale, only: MPR_UPSCALERS
  use mo_mpr_data_array_upscale, only: MPR_COORD_UPSCALERS
  use mo_mpr_parameters, only : MPR_PARAMETERS

  use mo_kind, only: i4

  implicit none

  private

  public :: reset
  ! --------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------
contains

  subroutine reset()
    !> deallocates all allocatable types and resets global vectors
    integer(i4) :: iItem

    ! also used in testing, so we need to check if allocated
    if (allocated(MPR_COORD_ALIAS)) deallocate(MPR_COORD_ALIAS)
    if (allocated(MPR_PARAMETERS%names)) deallocate(MPR_PARAMETERS%names)
    if (allocated(MPR_PARAMETERS%values)) deallocate(MPR_PARAMETERS%values)

    do iItem = 1, size(MPR_COORD_UPSCALERS)
      call MPR_COORD_UPSCALERS(iItem)%reset()
    end do

    do iItem = 1, size(MPR_UPSCALERS)
      call MPR_UPSCALERS(iItem)%reset()
    end do

    do iItem = 1, size(MPR_COORDINATES)
      call MPR_COORDINATES(iItem)%reset()
    end do

    do iItem = 1, size(MPR_DATA_ARRAYS)
      call MPR_DATA_ARRAYS(iItem)%reset()
    end do

  end subroutine reset

end module mo_mpr_reset

