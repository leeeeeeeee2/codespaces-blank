#include "flogging.h"
module mo_mpr_parameters

  use mo_kind, only: i4, dp
  use mo_constants, only: nodata_dp
  use mo_mpr_constants, only: defaultAlias, maxNameLength
  use mo_utils, only: eq
  use mo_orderpack, only: mrgrnk
  use flogging

  implicit none

  private

  public :: Parameters, mpr_add_parameters

  ! --------------------------------------------------------------------------------------
  type :: Parameters

    character(maxNameLength), dimension(:), allocatable :: names
    real(dp), dimension(:), allocatable :: values

  contains
    procedure, public :: does_contain
    procedure, public :: get_parameter

  end type Parameters

  type(Parameters), public :: MPR_PARAMETERS ! all coord_alias

  ! --------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------
contains

  function does_contain(self, name) result(is_contained)
    class(Parameters), intent(in) :: self
    character(*), intent(in) :: name
    logical :: is_contained

    is_contained = any(self%names == trim(name))
  end function does_contain

  function get_parameter(self, name) result(value)
    class(Parameters), intent(in) :: self
    character(*), intent(in) :: name
    real(dp) :: value

    integer(i4) :: i

    value = nodata_dp

    do i=1, size(self%names)
      if (self%names(i) == trim(name)) then
        value = self%values(i)
        exit
      end if
    end do

  end function get_parameter

  subroutine mpr_add_parameters(newNames, newValues)
    character(*), dimension(:), intent(in) :: newNames
    real(dp), dimension(:), intent(in) :: newValues

    character(maxNameLength), dimension(:), allocatable :: tempNames
    real(dp), dimension(:), allocatable :: tempValues

    integer(i4) :: n, counter

    ! check that both input arrays have the same size
    n = size(newNames)
    if (n /= size(newValues)) then
      log_error("(1X,A,A,I0,A,I0,A)") 'mpr_add_parameters: Provided ', n, ' parameter names and ', &
              size(newValues), ' parameter values.'
      stop 1
    end if
    ! add room for new values
    counter = 0
    ! check for state of allocation of the the MPR_PARAMETERS attributes
    if (allocated(MPR_PARAMETERS%names)) then
      ! they already contain data...
      ! copy the existing data in the temporary vectors and add room for new values
      counter = size(MPR_PARAMETERS%names)
    end if
    allocate(tempNames(size(newNames) + counter))
    allocate(tempValues(size(newNames) + counter))

    if (allocated(MPR_PARAMETERS%names)) then
      tempNames(1:counter) = MPR_PARAMETERS%names(:)
      tempValues(1:counter) = MPR_PARAMETERS%values(:)
      ! deallocate the existing values of the the MPR_PARAMETERS attributes
      deallocate(MPR_PARAMETERS%names)
      deallocate(MPR_PARAMETERS%values)
    end if
    tempNames(counter+1:size(tempNames)) = newNames
    tempValues(counter+1:size(tempValues)) = newValues

    ! remove duplicates in list to add when comparing with self and if applicable the existing values
    call remove_duplicates(tempNames, tempValues, MPR_PARAMETERS%names, MPR_PARAMETERS%values, &
            dropName=defaultAlias, dropValue=nodata_dp)
    deallocate(tempNames, tempValues)

  end subroutine mpr_add_parameters

  subroutine remove_duplicates(names, values, outValues, outJointValues, dropName, dropValue)
    !< removes duplicates in the array of names
    !< the items in the values array with the same indices of the duplicates are also removed
    !< additionally, items in both arrays are removed if an item of values equals dropValue

    character(*), dimension(:), intent(in) :: names
    real(dp), dimension(:), intent(in) :: values
    character(len(names)), dimension(:), intent(out), allocatable :: outValues
    real(dp), dimension(:), intent(out), allocatable :: outJointValues
    character(*), intent(in) :: dropName
    real(dp), intent(in) :: dropValue

    logical, dimension(size(names)) :: tempMask
    integer(i4), dimension(size(names)) :: tempIndices

    integer(i4) :: i, j, k

    if (size(names) == 0_i4) then
      allocate(outValues(0_i4))
      allocate(outJointValues(0_i4))
    else
      tempMask = .true.
      ! rank the names into the tempIndices array
      call mrgrnk(names, tempIndices)
      tempMask(1) = check_valid_parameter(trim(names(1)), values(1), trim(dropName), dropValue)
      do k = 2, size(tempIndices)
        i = tempIndices(k)
        j = tempIndices(k-1)
        ! check each item for dropValue and flag
        tempMask(i) = check_valid_parameter(trim(names(i)), values(i), trim(dropName), dropValue)
        ! check each item for similarity with its predecessor and flag
        if (tempMask(i) .and. names(i) == names(j)) then
          tempMask(i) = .false.
          log_warn(*) "you supplied a duplicate parameter name ('", trim(names(i)), &
                  "'), using value ", values(j), &
                  " and dropping ", values(i), "."
        end if
      end do
      outValues = pack(names, mask=tempMask)
      outJointValues = pack(values, mask=tempMask)
    end if

  end subroutine remove_duplicates
  
  function check_valid_parameter(name, value, dropName, dropValue) result(boolean)
    character(*), intent(in) :: name
    real(dp), intent(in) :: value
    character(*), intent(in) :: dropName
    real(dp), intent(in) :: dropValue
    logical :: boolean

    boolean = .true.
    if (name == dropName) then
      boolean = .false.
    else if (eq(value, dropValue)) then
      boolean = .false.
      log_warn(*) "you supplied a parameter name ('", name, "'), but no value."
    end if
    
  end function check_valid_parameter

end module mo_mpr_parameters

