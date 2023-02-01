#include "flogging.h"
module mo_mpr_upscale_func
  use mo_kind, only: dp, i4
  use mo_utils, only: eq
  !use ieee_arithmetic, only : ieee_is_nan
  use mo_constants, only: nodata_dp
  use flogging

  implicit none

  private

  public :: wrap_upscale, wrap_weighted_upscale, get_upscale_func, upscale_func_alias

  abstract interface
    real(dp) function upscale_func_alias(array, weights, p)
      ! import the double precision kind specification
      import dp
      !> the array of subgrid values (no missing values)
      real(dp), dimension(:), intent(in) :: array
      !> the array of weights (same shape as array)
      real(dp), dimension(:), intent(in), optional :: weights
      !> an optional parameter passed to the function (power mean)
      real(dp), intent(in), optional :: p
    end function upscale_func_alias
  end interface

contains

  ! ----------------------------------------------------------------------------------------
  subroutine get_upscale_func(upscaleOperatorName, upscaleFunc, p, doNeedWeights)
    !< use a string denoting the upscale func and return a pointer to the correct function and also the p-norm

    !> this is the input string
    character(*), intent(in) :: upscaleOperatorName
    !> pointer to an upscale_func_alias procedure
    procedure(upscale_func_alias), pointer, intent(out) :: upscaleFunc
    !> p-norm parameter
    real(dp), intent(out) :: p
    !> flag whether to point to functions using weights
    logical, intent(in) :: doNeedWeights

    ! correctly initialize
    upscaleFunc => null()
    ! check for first character
    select case(upscaleOperatorName(:1))
      ! if it is a letter, then assume a string
      case('a':'z')
        p = nodata_dp
        select case(trim(upscaleOperatorName))
          case('sum')
            upscaleFunc => sum_func
          case('min')
            upscaleFunc => min_func
          case('max')
            upscaleFunc => max_func
          case('var')
            if (doNeedWeights) then
              upscaleFunc => weighted_var_func
            else
              upscaleFunc => var_func
            end if
          case('std')
            if (doNeedWeights) then
              upscaleFunc => weighted_std_func
            else
              upscaleFunc => std_func
            end if
          case('laf')
            upscaleFunc => laf_func
          case default
            log_error(*) 'Upscaling operator "', trim(upscaleOperatorName), '" is not valid.', &
                    ' Use one of "sum", "min", "max", "var", "std", "laf" or a floating point number'
            stop 1
        end select
      ! if it is a number, then assume generalized mean
      case('0':'9', '-')
        ! assume a real number
        ! TODO: How to capture possible errors?
        read(upscaleOperatorName, *) p
        if (eq(p, 0.0_dp)) then
          ! geometric mean special case
            if (doNeedWeights) then
              upscaleFunc => weighted_geometric_mean
            else
              upscaleFunc => geometric_mean
            end if
        else
          ! all other cases
            if (doNeedWeights) then
              upscaleFunc => weighted_generalized_mean
            else
              upscaleFunc => generalized_mean
            end if
        end if
      ! all other cases are not allowed
      case default
        log_error(*) 'Upscaling operator "', trim(upscaleOperatorName), '" is not valid.', &
                ' Use one of "sum", "min", "max", "var", "std", "laf" or a floating point number'
        stop 1
    end select

  end subroutine get_upscale_func

  subroutine wrap_weighted_upscale(func, p, weightsIn, sliceIn, maskIn, nanWeightsIn, valueOut, maskOut, nanWeightsOut)
    !< wrapper function for applying the upscaling function
    !< it handles possible existing missing values and updates the weights accordingly
    procedure(upscale_func_alias), intent(in), pointer :: func
    real(dp), intent(in) :: p
    real(dp), dimension(:), intent(in)  :: weightsIn
    real(dp), dimension(:), intent(in)  :: sliceIn
    logical, dimension(:), intent(in)  :: maskIn
    real(dp), dimension(:), intent(in), optional  :: nanWeightsIn
    real(dp), intent(out) :: valueOut
    logical, intent(out) :: maskOut
    real(dp), intent(out), optional :: nanWeightsOut

    real(dp), dimension(size(weightsIn)) :: weightsAdapted

    if (all(.not. maskIn)) then
      ! if all values are masked missing, then the result is nan and also masked missing
      valueOut = nodata_dp
      maskOut = .false.
      if (present(nanWeightsOut)) nanWeightsOut = 0.0_dp
    else
      if (present(nanWeightsOut) .and. present(nanWeightsIn)) then
        ! the weights based on the previous upscaling are combined with the weights based on different cell widths
        weightsAdapted = nanWeightsIn * weightsIn
        ! the resulting nan_weight is the sum of all valid cells
        nanWeightsOut = sum(weightsAdapted, mask=maskIn)
        ! they are rescaled so they sum to 1
        weightsAdapted = weightsAdapted / nanWeightsOut
      else
        weightsAdapted = weightsIn / sum(weightsIn, mask=maskIn)
      end if
      maskOut = .true.
      ! the non-missing values are send to func together with their weights
      valueOut = func(pack(sliceIn, mask=maskIn), pack(weightsAdapted, mask=maskIn), p)
    end if

  end subroutine wrap_weighted_upscale

  subroutine wrap_upscale(func, p, sliceIn, maskIn, valueOut, maskOut)
    !< lightweight wrapper function for applying the upscaling function (in case of no weights necessary)
    procedure(upscale_func_alias), intent(in), pointer :: func
    real(dp), intent(in) :: p
    real(dp), dimension(:), intent(in)  :: sliceIn
    logical, dimension(:), intent(in)  :: maskIn
    real(dp), intent(out) :: valueOut
    logical, intent(out) :: maskOut

    if (all(.not. maskIn)) then
      ! if all values are masked missing, then the result is nan and also masked missing
      valueOut = nodata_dp
      maskOut = .false.
    else
      maskOut = .true.
      valueOut = func(pack(sliceIn, mask=maskIn), p=p)
    end if

  end subroutine wrap_upscale

  real(dp) function generalized_mean(array, weights, p)
    !< calculate the generalized mean of an array
    real(dp), dimension(:), intent(in) :: array
    real(dp), dimension(:), intent(in), optional :: weights
    real(dp), intent(in), optional :: p

    generalized_mean = (sum(array**p)/size(array))**(1.0_dp / p)

  end function generalized_mean

  real(dp) function geometric_mean(array, weights, p)
    !< calculate the geometric mean of an array
    real(dp), dimension(:), intent(in) :: array
    real(dp), dimension(:), intent(in), optional :: weights
    real(dp), intent(in), optional :: p

    geometric_mean = product(array) ** (1.0_dp / size(array))

  end function geometric_mean

  real(dp) function sum_func(array, weights, p)
    !< calculate the sum of an array
    real(dp), dimension(:), intent(in) :: array
    real(dp), dimension(:), intent(in), optional :: weights
    real(dp), intent(in), optional :: p

    sum_func = sum(array)

  end function sum_func

  real(dp) function min_func(array, weights, p)
    !< get the minimum of an array
    real(dp), dimension(:), intent(in) :: array
    real(dp), dimension(:), intent(in), optional :: weights
    real(dp), intent(in), optional :: p

    min_func = minval(array)

  end function min_func

  real(dp) function max_func(array, weights, p)
    !< get the maximum of an array
    real(dp), dimension(:), intent(in) :: array
    real(dp), dimension(:), intent(in), optional :: weights
    real(dp), intent(in), optional :: p

    max_func = maxval(array)

  end function max_func

  real(dp) function var_func(array, weights, p)
    !< calculate the weighted variance of an array
    real(dp), dimension(:), intent(in) :: array
    real(dp), dimension(:), intent(in), optional :: weights
    real(dp), intent(in), optional :: p

    real(dp) :: avg
    integer(i4) :: n

    n = size(array)
    avg = sum(array) / n
    var_func = sum((array - avg) ** 2.0_dp) / n

  end function var_func

  real(dp) function std_func(array, weights, p)
    !< calculate the weighted standard deviation of an array
    real(dp), dimension(:), intent(in) :: array
    real(dp), dimension(:), intent(in), optional :: weights
    real(dp), intent(in), optional :: p

    real(dp) :: avg
    integer(i4) :: n

    n = size(array)
    avg = sum(array) / n
    std_func = (sum((array - avg) ** 2.0_dp) / n) ** 0.5_dp

  end function std_func

  real(dp) function weighted_generalized_mean(array, weights, p)
    !< calculate the generalized mean of an array
    real(dp), dimension(:), intent(in) :: array
    real(dp), dimension(:), intent(in), optional :: weights
    real(dp), intent(in), optional :: p

    weighted_generalized_mean = sum(weights * array**p)**(1.0_dp / p)

  end function weighted_generalized_mean

  real(dp) function weighted_geometric_mean(array, weights, p)
    !< calculate the geometric mean of an array
    real(dp), dimension(:), intent(in) :: array
    real(dp), dimension(:), intent(in), optional :: weights
    real(dp), intent(in), optional :: p

    weighted_geometric_mean = product(array ** weights)

  end function weighted_geometric_mean

  real(dp) function weighted_var_func(array, weights, p)
    !< calculate the weighted variance of an array
    real(dp), dimension(:), intent(in) :: array
    real(dp), dimension(:), intent(in), optional :: weights
    real(dp), intent(in), optional :: p

    real(dp) :: avg

    avg = sum(array * weights)
    weighted_var_func = sum(weights *((array - avg) ** 2.0_dp))

  end function weighted_var_func

  real(dp) function weighted_std_func(array, weights, p)
    !< calculate the weighted standard deviation of an array
    real(dp), dimension(:), intent(in) :: array
    real(dp), dimension(:), intent(in), optional :: weights
    real(dp), intent(in), optional :: p

    real(dp) :: avg

    avg = sum(array * weights)
    weighted_std_func = (sum(weights *((array - avg) ** 2.0_dp))) ** 0.5_dp

  end function weighted_std_func

  real(dp) function laf_func(array, weights, p)
    real(dp), dimension(:), intent(in) :: array
    real(dp), dimension(:), intent(in), optional :: weights
    real(dp), intent(in), optional :: p

    integer(i4), dimension(1) :: ind

    ind = maxloc(weights)

    laf_func = array(ind(1))

  end function laf_func

end module mo_mpr_upscale_func
