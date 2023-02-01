#include "flogging.h"
module mo_mpr_transfer_func
  use mo_kind, only : dp, i8, i4
  use mo_mpr_input_field_container, only: InputFieldContainer
  use mo_constants, only: nodata_dp
  use mo_mpr_constants, only: defaultAlias, maxStringLength
  use mo_orderpack, only: sort_index
  use mo_utils, only: eq
  use flogging

  implicit none

  private

  public :: transfer_function_2
  public :: transfer_function_1

  type, public :: TransferFunctionTable
    integer(i4), dimension(3) :: indices = [1_i4, 1_i4, 2_i4]
    character(maxStringLength), dimension(3) :: names = [character(maxStringLength) :: &
            "identity", defaultAlias, "empirical_cdf"]
  contains
    private
    procedure, public :: get_index

  end type TransferFunctionTable

  type(TransferFunctionTable), public :: mprTransferFunctionTable

  abstract interface
    function transfer_func_alias(x, param) result(func_result)
      ! import the double precision kind specification and custom type
      import dp, InputFieldContainer
      !> an array containing the predictor variables (access values through `data_p` property)
      type(InputFieldContainer), intent(in) :: x(:)
      !> an array containing the TF parameters
      real(dp), intent(in) :: param(:)
      !> the resulting TF result
      real(dp), allocatable :: func_result(:)
    end function transfer_func_alias
  end interface

contains
  function get_index(self, name) result(index)
    class(TransferFunctionTable), intent(in):: self
    character(*), intent(in):: name
    integer(i4) :: index, iName

    do iName=1, size(self%names)
      if (trim(self%names(iName)) == trim(name)) then
        index = self%indices(iName)
        return
      end if
    end do

    log_warn(*) 'The transfer function name "', trim(name) , '" is not registered.'
    index = 0_i4

  end function get_index

  function transfer_function_1(x, param) result(func_result)
    real(dp), dimension(:), allocatable :: func_result
    type(InputFieldContainer), dimension(:), intent(in) :: x
    real(dp), dimension(:), intent(in) :: param
    integer(i8) :: n

    n = size(x(1)%data_p)
    allocate(func_result(n))
    func_result = x(1)%data_p
  end function transfer_function_1

  function transfer_function_2(x, param) result(func_result)
    real(dp), dimension(:), allocatable :: func_result
    type(InputFieldContainer), dimension(:), intent(in) :: x
    real(dp), dimension(:), intent(in) :: param

    integer(i4), dimension(:), allocatable :: valueSortedIndex
    integer(i8) :: n, i, iSort, iSortAfter

    n = size(x(1)%data_p, kind=i8)
    allocate(func_result(n), valueSortedIndex(n))

    ! get sorted data and sorted indexes to remap later
    valueSortedIndex = sort_index(x(1)%data_p)

    ! empirical distribution of values = cumulated number points with values that are <= the value at this point
    !
    !       sorted data                     emp. CDF
    ! 9 |             x x       7/8 |             x x
    !   |                           |
    ! 8 |           x           5/8 |           x
    !   |                           |
    ! 5 |     x x x             4/8 |     x x x
    !   |                           |
    ! 2 |  x                    1/8 |  x
    !   |__________________         |__________________
    !
    ! EXAMPLE
    ! in      = [  7, 20, 31, 31, 12, 31, 42 ]
    ! sorted  = [  7, 12, 20, 31, 31, 31, 42 ]
    ! index   = [  1,  5,  2,  3,  4,  6,  7 ]
    ! temp    = [  1,  2,  3,  6,  6,  6,  7 ]
    ! out     = [  1,  3,  6,  6,  2,  6,  7 ] / (len(out) + 1 )

    ! highest value value = highest rank or No. of data points / (data points + 1)
    func_result(valueSortedIndex(n)) = real(n, dp) / real(n + 1_i4, dp)

    ! backward loop to check if the preceding data point has the same value value
    do i = n - 1, 1, -1
      iSort=valueSortedIndex(i)
      iSortAfter=valueSortedIndex(i+1)
      if (eq(x(1)%data_p(iSort), x(1)%data_p(iSortAfter))) then
        ! if yes: assign the same probabitity
        func_result(iSort) = func_result(iSortAfter)
      else
        ! if not: assign rank / (data points + 1)
        func_result(iSort) = real(i, dp) / real(n + 1_i4, dp)
      end if
    end do

    deallocate(valueSortedIndex)

  end function transfer_function_2
end module mo_mpr_transfer_func
