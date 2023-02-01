#include "flogging.h"
module mo_mpr_data_array_upscale

  use mo_kind, only : i4, i8, dp
  use mo_mpr_constants, only : maxNoDataArrays, defaultMapMethod, scripDstAddressName,&
          scripSrcAddressName, scripWeightsName, scripLinksCoord, scripWeightsCoord, defaultCoordCount, &
          maxTolerance
  use mo_constants, only : nodata_dp, nodata_i4, eps_dp
  use mo_mpr_coordinate, only : Coordinate, CoordinatePointer, get_index_in_coordinate, &
          create_upscaler_name, create_upscaler_alias, combine_coordinates, create_integral_coordinate, &
          split_coordinate, check_within_bounds
  use mo_mpr_util_types, only : MprBaseType
  use mo_mpr_utils, only : get_n_initialized
  use flogging
  use mo_append, only : append
  use mo_mpr_constants, only : maxNameLength, maxStringLength, defaultAlias
  use mo_utils, only: ne, eq, ge, le, unpack_chunkwise
  use mo_mpr_upscale_func, only : wrap_upscale, wrap_weighted_upscale, get_upscale_func, upscale_func_alias
  use mo_netcdf, only : NcDataset, NcVariable, NcDimension
  use mo_poly, only: inpoly, orientpoly, mod_pole, mod_shift

  implicit none

  private

  public :: add_upscaler
  public :: add_coord_upscaler, get_index_in_coord_upscaler
  public :: get_target_ids, get_coordUpscaler, compute_index_helpers, permutate

  ! --------------------------------------------------------------------------------------
  type, extends(MprBaseType), public :: UpscaleHelper
    private
    !> flag whether to execute upscaling
    logical, public :: doUpscale
    !> the source coordinates (they will be altered by the init function according to reordering and combination
    type(CoordinatePointer), dimension(:), allocatable :: sourceCoords
    !> the target coordinates
    type(CoordinatePointer), dimension(:), allocatable :: targetCoords
    !> flag whether to execute reordering of coordinates
    logical :: doReorder
    !> reorder indices based on original order
    integer(i4), dimension(:), allocatable :: reorderIndices
    !> original shape of array (size of coordinates)
    integer(i8), dimension(:), allocatable :: originalShape
    !> flag whether to execute reordering of coordinates
    logical :: doSplitTargetCoord
    !> flag whether combining of source coords happened
    logical :: didCombineSourceCoords
    !> indices of (temporary) target coordinate to split into elementary coordinates
    integer(i4), dimension(:), allocatable :: splitTargetIndices
    !> indices of (original) sources coordinate combined into 2d coordinate
    integer(i4), dimension(:), allocatable :: combineSourceIndices
    !> newShapes stores the intermediate shapes during upscaling
    integer(i8), dimension(:,:), allocatable :: newShapes
    !> index of CoordUpscalers to use for upscaling
    integer(i4), dimension(:), allocatable :: upscalersIndex
    !> whether to upscale that coordinate
    logical, dimension(:), allocatable :: doCoordUpscale
    !> strings of upscale operators
    character(maxNameLength), dimension(:), allocatable :: upscaleOperatorNames


  contains
    private
    procedure, public :: is_finalized => is_finalized_UpscaleHelper
    procedure, public :: reset =>  reset_UpscaleHelper
    procedure, public :: execute => execute_UpscaleHelper
    procedure, public :: get_stats => get_stats_UpscaleHelper
    procedure :: prepare_coords
    procedure :: combine_coords
    procedure :: check_order_coords
    procedure :: check_split_target_dim
    procedure :: split_target_coordinate
    procedure :: trim_upscale_operators
    procedure :: init_remaining_target_dim

  end type UpscaleHelper

  interface UpscaleHelper
    procedure newUpscaleHelper
  end interface UpscaleHelper

  type(UpscaleHelper), dimension(maxNoDataArrays), public, target :: MPR_UPSCALERS

  type, extends(MprBaseType), public :: CoordUpscaler
    private
    character(maxNameLength), public :: alias
    integer(i8), dimension(:), allocatable :: subcells
    integer(i8), dimension(:, :), allocatable :: ids
    real(dp), dimension(:, :), allocatable :: weights
    character(maxStringLength), public :: mapMethod

    logical :: doNeedWeights

  contains
    private
    procedure :: compute_weights_1d
    procedure :: compute_weights_1d_multiple
    procedure :: compute_weights_1d_
    procedure :: compute_weights_2d
    procedure :: compute_weights_poly
    procedure :: calculate_weights
    procedure :: combine_weights
    procedure :: read_weights
    procedure :: write_weights
    procedure :: from_other
    procedure :: init_values
    procedure :: check_step_multiples
    procedure, public :: get_stats
    procedure, public :: is_finalized => is_finalized_CoordUpscaler
    procedure, public :: reset => reset_CoordUpscaler
    procedure, public :: execute => execute_CoordUpscaler
    procedure, public :: get_n_target_cells
    procedure, public :: get_max_subcells
    procedure, public :: get_ids
    procedure, public :: get_subcells
    procedure, public :: get_weights

  end type CoordUpscaler

  interface CoordUpscaler
    procedure newCoordUpscaler
  end interface CoordUpscaler

  type(CoordUpscaler), dimension(maxNoDataArrays), public, target :: MPR_COORD_UPSCALERS

  type, public :: CoordUpscalerPointer
    ! this is a pointer to a Coordinate
    type(CoordUpscaler), pointer :: coord_p => null()
  contains
    private

  end type CoordUpscalerPointer

  abstract interface
     function get_source_index_alias(targetBound, vals, cachedIndex, step) result(index)
      !< find index of lower vals where a comparison with targetBound is not met anymore
      import dp, i8
      !> bound value to compare values against
      real(dp), intent(in) :: targetBound
      !> values to be compared (ordered descendingly, unequal step size)
      real(dp), dimension(:), intent(in) :: vals
      !> optional starting point for index search
      integer(i8), intent(in) :: cachedIndex
      !> step size of vals, optional and in this case never used
      real(dp), intent(in), optional :: step
      !> return value: index of vals
      integer(i8) :: index
    end function get_source_index_alias
  end interface

  ! --------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------
contains

  function newUpscaleHelper(name, id, doUpscale, targetCoords, sourceCoords, upscaleOperatorNames)
    character(*), intent(in) :: name
    integer(i4), intent(in) :: id
    logical, intent(in), optional :: doUpscale
    type(CoordinatePointer), dimension(:), intent(in), optional :: targetCoords
    type(CoordinatePointer), dimension(:), intent(in), optional :: sourceCoords
    character(maxNameLength), dimension(:), intent(in), optional :: upscaleOperatorNames
    type(UpscaleHelper) :: newUpscaleHelper

    newUpscaleHelper%name = trim(name)
    newUpscaleHelper%id = id
    newUpscaleHelper%is_initialized = .true.
    newUpscaleHelper%doReorder = .false.
    allocate(newUpscaleHelper%reorderIndices(0))
    newUpscaleHelper%doSplitTargetCoord = .false.
    allocate(newUpscaleHelper%splitTargetIndices(0))
    newUpscaleHelper%didCombineSourceCoords = .false.
    allocate(newUpscaleHelper%combineSourceIndices(0))

    ! check for main flag
    if (present(doUpscale)) then
      newUpscaleHelper%doUpscale = doUpscale
    else
      newUpscaleHelper%doUpscale = .true.
    end if

    ! if main flag is okay, check for existance of sourceCoords, targetCoords
    if (newUpscaleHelper%doUpscale) then
      if (.not. present(targetCoords)) then
        log_error(*) 'Target coordinates not provided for UpscaleHelper', trim(newUpscaleHelper%name)
        stop 1
      end if

      if (.not. present(sourceCoords)) then
        log_error(*) 'Source coordinates not provided for UpscaleHelper', trim(newUpscaleHelper%name)
        stop 1
      end if

      if (present(upscaleOperatorNames)) then
        newUpscaleHelper%upscaleOperatorNames = upscaleOperatorNames
      else
        log_error(*) 'upscaleOperatorNames not provided for UpscaleHelper', trim(newUpscaleHelper%name)
        stop 1
      end if

      ! check which operations might be performed and alter coordinates
      call newUpscaleHelper%prepare_coords(sourceCoords, targetCoords)
    end if

  end function newUpscaleHelper

  function is_finalized_UpscaleHelper(self) result(is_finalized)
    class(UpscaleHelper), intent(in) :: self
    logical :: is_finalized

    if (self%doUpscale) then
      is_finalized = allocated(self%targetCoords)
    else
      is_finalized = .true.
    end if
  end function is_finalized_UpscaleHelper

  subroutine  reset_UpscaleHelper(self)
    class(UpscaleHelper), intent(inout) :: self

    if (self%is_initialized) then
      self%name = ''
      self%id=0
      self%is_initialized = .false.
      if (allocated(self%sourceCoords)) deallocate(self%sourceCoords)
      if (allocated(self%targetCoords)) deallocate(self%targetCoords)
    end if

  end subroutine  reset_UpscaleHelper

  subroutine prepare_coords(self, sourceCoords, targetCoords)
    class(UpscaleHelper) :: self
    type(CoordinatePointer), dimension(:), intent(in) :: sourceCoords
    type(CoordinatePointer), dimension(:), intent(in) :: targetCoords

    ! check the order of indices
    call self%check_order_coords(sourceCoords, targetCoords)

    ! we need to check if any found indices are occurring twice - this means that a xy source coordinate
    ! needs to be upscaled into x and y coordinates
    ! -> we temporarily replace the x and y target by xy and schedule a split of coordinates for after upscaling
    call self%check_split_target_dim(sourceCoords, targetCoords)

    ! init the remaining targetCoord vectors if this has not yet happened
    call self%init_remaining_target_dim()

    if (self%doSplitTargetCoord) then
      ! trim upscaleOperatorNames with self%splitTargetIndices, because we know that targetCoords was trimmed
      ! after reading config
      call self%trim_upscale_operators()
    end if

    ! sanity check
    if (size(self%upscaleOperatorNames) /= size(self%targetCoords)) then
      log_error(*) "Error initializing UpscaleHelper: '", trim(self%name), &
              "' received ", size(self%upscaleOperatorNames), &
              " upscale operators and ", size(self%targetCoords), " coordinates. Number must be equal."
      call self%get_stats()
      stop 1
    end if
    ! combine coords based on same upscale operator
    ! e.g. upscalers (for target coords 1, 2, 3, 4) with upscale_ops 1,1,1,1   -> new upscalers 5
    ! e.g. upscalers (for target coords 1, 2, 3, 4) with upscale_ops 1,-1,-1,1 -> new upscalers 1,5,4
    ! e.g. upscalers (for target coords 1, 2, 3, 4) with upscale_ops -1,-1,1,1 -> new upscalers 5,6
    call self%combine_coords()


  end subroutine prepare_coords

  subroutine combine_coords(self)
    class(UpscaleHelper) :: self

    character(maxNameLength), dimension(:), allocatable :: tempUpscaleOperatorNames
    integer(i4), dimension(:), allocatable :: tempUpscalersIndex
    integer(i4) :: iCoord, nCoord, id, iCoordCurrent, slots, jCoord
    logical :: doCombineUpscalers
    integer(i8), dimension(:), allocatable :: sourceShape, targetShape

    ! init the values for the properties
    nCoord = size(self%targetCoords)
    allocate(self%doCoordUpscale(nCoord))
    ! set it first to maximum possible value 
    allocate(self%upscalersIndex(nCoord))
    self%upscalersIndex = nodata_i4
    self%doCoordUpscale = .false.

    allocate(sourceShape(nCoord))
    allocate(targetShape(nCoord))
    allocate(tempUpscaleOperatorNames(nCoord))

    ! counter for all upscalers needed
    slots = 0_i4
    doCombineUpscalers = .true.
    ! init counter of current iCoord
    iCoordCurrent = 1
    do iCoord = 1, nCoord
      ! advance counter only, if we do not have a following upscaling operator with the same operator
      if (doCombineUpscalers) iCoordCurrent = iCoord
      ! check if upscaling is necessary because unequal coordinates
      if (self%sourceCoords(iCoord)%coord_p /= self%targetCoords(iCoord)%coord_p) then
        self%doCoordUpscale(iCoord) = .true.
        ! check if we are the last coordinate
        if (iCoord + 1 > nCoord) then
          doCombineUpscalers = .true.
        ! look ahead if next coordinate also requires upscaling and has similar upscaling operator
        ! if not then create
        else if (self%sourceCoords(iCoord+1)%coord_p /= self%targetCoords(iCoord+1)%coord_p .and. &
          self%upscaleOperatorNames(iCoord) == self%upscaleOperatorNames(iCoord+1)) then
          doCombineUpscalers = .false.
        else
          doCombineUpscalers = .true.
        end if
      end if
      if (self%doCoordUpscale(iCoord) .and. doCombineUpscalers) then
        ! increase counter
        slots = slots + 1_i4
        ! create upscaler
        call get_coordUpscaler(self%sourceCoords(iCoordCurrent:iCoord), self%targetCoords(iCoordCurrent:iCoord), id)
        log_trace(*) 'combine_coords: created upscaler with id', id, trim(MPR_COORD_UPSCALERS(id)%name)
        self%upscalersIndex(slots) = id
        tempUpscaleOperatorNames(slots) = self%upscaleOperatorNames(iCoord)
        sourceShape(slots) = product([(self%sourceCoords(jCoord)%coord_p%count, &
                jCoord=iCoordCurrent, iCoord)])
        targetShape(slots) = product([(self%targetCoords(jCoord)%coord_p%count, &
                jCoord=iCoordCurrent, iCoord)])
        log_trace(*) 'combine_coords: using jCoords from ', iCoordCurrent, ' to ', iCoord, &
                ' for the sourceShape with the product of these counts:', &
                [(self%sourceCoords(jCoord)%coord_p%count, jCoord=iCoordCurrent, iCoord)], ': ',&
                sourceShape(slots)
        log_trace(*) 'combine_coords: using jCoords from ', iCoordCurrent, ' to ', iCoord, &
                ' for the targetShape with the product of these counts:', &
                [(self%targetCoords(jCoord)%coord_p%count, jCoord=iCoordCurrent, iCoord)], ': ',&
                targetShape(slots)
        log_debug(*) 'combine_coords: combined upscaler ', trim(MPR_COORD_UPSCALERS(id)%name), &
              ' which is now used for target coordinates at positions ', iCoordCurrent, ':', iCoord
      else if (.not. self%doCoordUpscale(iCoord)) then
        ! increase counter
        slots = slots + 1_i4
        self%upscalersIndex(slots) = nodata_i4
        tempUpscaleOperatorNames(slots) = defaultAlias
        sourceShape(slots) = self%sourceCoords(iCoord)%coord_p%count
        targetShape(slots) = self%targetCoords(iCoord)%coord_p%count
        log_debug(*) 'combine_coords: no upscaling necessary for iCoord ', iCoord, &
                ' with name ', trim(self%targetCoords(iCoord)%coord_p%name)
      end if
    end do

    ! reduce array self%upscalersIndex size to needed length
    allocate(tempUpscalersIndex(slots))
    tempUpscalersIndex = self%upscalersIndex(1:slots)
    call move_alloc(tempUpscalersIndex, self%upscalersIndex)
    ! reduce array self%upscaleOperatorsNames size to needed length
    deallocate(self%upscaleOperatorNames)
    allocate(self%upscaleOperatorNames(slots))
    self%upscaleOperatorNames = tempUpscaleOperatorNames(1:slots)
    deallocate(tempUpscaleOperatorNames)

    ! what newShapes is: each row contains an intermediate shape of the array for upscaling purposes
    ! the first row contains the final state, the second the previous and the last the current (initial) shape
    ! thus there are slots + 1 rows to cover all states
    allocate(self%newShapes(slots, slots + 1))
    ! calculate intermediate shapes of the newArray
    do iCoord = 1, slots
      self%newShapes(slots+1-iCoord, :) = sourceShape(iCoord)
      self%newShapes(slots+1-iCoord, iCoord+1:slots+1) = targetShape(iCoord)
    end do
    log_debug(*) 'combine_coords: the newShapes array containing all array shapes for each ', &
            'intermediate upscaling steps is (shape: ', shape(self%newShapes), '):', self%newShapes

    if (all(.not. self%doCoordUpscale)) then
      self%doUpscale = .false.
    end if

  end subroutine combine_coords


  subroutine check_order_coords(self, sourceCoords, targetCoords)
    class(UpscaleHelper) :: self
    type(CoordinatePointer), dimension(:), intent(in) :: sourceCoords
    type(CoordinatePointer), dimension(:), intent(in) :: targetCoords

    integer(i4) :: iCoord, jCoord, prevCoordIndex, id
    integer(i4), dimension(:), allocatable :: foundIndices

    ! check the order of indices
    allocate(self%sourceCoords(size(targetCoords)))
    do iCoord = 1, size(targetCoords)
      ! this returns the indices in the right order, if multiple coordinates need to be combined
      foundIndices = get_index_in_coordinate(&
              coordinateName=targetCoords(iCoord)%coord_p%name, &
              coordinates=sourceCoords, &
              checkAliasesArg=.true.)
      log_debug(*) 'check_order_coords: target coordinate name ', trim(targetCoords(iCoord)%coord_p%name), &
        ' got indices ', foundIndices, ' when checking coordinate (', &
        [(sourceCoords(jCoord)%coord_p%name, jCoord=1, size(sourceCoords))], ')'

      if (size(foundIndices) > 1_i4) then
        ! the source coordinates need to be combined to a new coordinate
        call combine_coordinates(sourceCoords(foundIndices), id)
        ! reference it
        call self%sourceCoords(iCoord)%set_coordinate_pointer(id)
        ! set flags
        self%didCombineSourceCoords = .true.
        self%combineSourceIndices = foundIndices
      else if (foundIndices(1) > size(sourceCoords)) then
        ! broadcasting is necessary
        ! create a new coordinate
        call create_integral_coordinate(targetCoords(iCoord)%coord_p, id)
        ! reference it
        call self%sourceCoords(iCoord)%set_coordinate_pointer(id)
      else
        self%sourceCoords(iCoord) = sourceCoords(foundIndices(1))
      end if
      call append(self%reorderIndices, foundIndices)
    end do

    ! now we need to check if the number of targetCoords >= sourceCoords
    do iCoord=1, size(sourceCoords)
      ! check if there is a source coordinate, that is not appearing in the reordered indices
      if (count(self%reorderIndices == iCoord) == 0_i4) then
        log_error(*) 'The coordinate ', trim(sourceCoords(iCoord)%coord_p%name), &
                ' needs to be upscaled to a target coordinate, but you did not provide any.'
        call self%get_stats()
        stop 1
      end if
    end do

    ! we now need to check, if the reorderIndices are ordered
    prevCoordIndex = 0_i4
    do iCoord = 1, size(self%reorderIndices)
      ! check if this is a coordinate that is not arbitrarily input ("broadcasting")
      if (self%reorderIndices(iCoord) <= size(sourceCoords)) then
        if (self%reorderIndices(iCoord) < prevCoordIndex) then
          ! this is very expensive, issue a warning, so the user might reorder their input data
          log_warn(*) 'The coordinate ', trim(targetCoords(iCoord)%coord_p%name), ' is not in the right order.', &
              ' Check, if a reordering of the input arrays is feasible to avoiding expensive transpose operations.'
          ! we need to reorder
          self%doReorder = .true.
        end if
        prevCoordIndex = self%reorderIndices(iCoord)
      end if
    end do

    ! prepare the remaining values for reordering (storing the original sizes)
    if (self%doReorder) then
      allocate(self%originalShape(size(sourceCoords)))
      do iCoord = 1, size(sourceCoords)
        self%originalShape(iCoord) = sourceCoords(iCoord)%coord_p%count
      end do
    end if

  end subroutine check_order_coords


  subroutine check_split_target_dim(self, sourceCoords, targetCoords)
    class(UpscaleHelper) :: self
    type(CoordinatePointer), dimension(:), intent(in) :: sourceCoords
    type(CoordinatePointer), dimension(:), intent(in) :: targetCoords

    integer(i4) :: id
    integer(i8) :: iCoord, jCoord, indexCounter
    logical, dimension(:), allocatable :: mask
    type(CoordinatePointer), dimension(:), allocatable :: tempCoordPointers

    ! we need to check if any found indices are occurring twice - this means that a xy source coordinate
    ! needs to be upscaled into x and y coordinates
    ! -> we temporarily replace the x and y target by xy and schedule a split of coordinates for after upscaling
    do iCoord = 1, size(self%reorderIndices)
      ! only use indices actually present in source (ignore broadcasting coordinates)
      if (self%reorderIndices(iCoord) > size(sourceCoords)) cycle
      ! check if an index occurs more than once in the remaining values
      mask = (self%reorderIndices(iCoord:) == self%reorderIndices(iCoord))
      if (count(mask, kind=i8) > 1_i8) then
        ! check if we have already found multiple indices before
        if (self%doSplitTargetCoord) then
          log_error(*) 'We do not yet support splitting of multiple target coordinates.', &
              'This is needed for UpscaleHelper ', trim(self%name)
          call self%get_stats()
          stop 1
        end if
        self%doSplitTargetCoord = .true.
        ! init the array with the proper length
        self%splitTargetIndices = pack(self%reorderIndices(iCoord:), mask)
        self%splitTargetIndices(1) = int(iCoord, kind=i4)
        indexCounter = 2_i8
        do jCoord=iCoord+1_i8, size(mask, kind=i8)
          ! get the indices of all the splitTargets
          if (mask(jCoord)) then
            self%splitTargetIndices(indexCounter) = int(jCoord, kind=i4)
            ! check if they are neighboring Indices, we do not yet support post-upscale reordering of coordinates
            if (self%splitTargetIndices(indexCounter) - self%splitTargetIndices(indexCounter-1) > 1_i4) then
              log_error(*) 'We do not yet support post-upscale reordering of coordinates.', &
                      'This is needed for UpscaleHelper ', trim(self%name)
              call self%get_stats()
              stop 1
            end if
            ! advance counter
            indexCounter = indexCounter + 1_i4
          end if
        end do
        ! now alter the target_coordinates
        call combine_coordinates(targetCoords(self%splitTargetIndices), id)
      end if
    end do

    if (self%doSplitTargetCoord) then
      ! allocate the proper target coords
      allocate(self%targetCoords(size(targetCoords) -1))
      ! also modify the self%sourceCoords as they still have a duplicate entry
      allocate(tempCoordPointers(size(targetCoords) -1))
      ! and set the correct pointers
      indexCounter = 0_i4
      do iCoord=1, size(self%targetCoords)
        tempCoordPointers(iCoord) = self%sourceCoords(iCoord + indexCounter)
        if (self%splitTargetIndices(1) == iCoord) then
          ! let the coordinate point to our new coordinate
          call self%targetCoords(iCoord)%set_coordinate_pointer(id)
          indexCounter = indexCounter + 1_i4
        else
          self%targetCoords(iCoord) = targetCoords(iCoord + indexCounter)
        end if
      end do
      call move_alloc(tempCoordPointers, self%sourceCoords)
    else
      self%targetCoords = targetCoords
    end if

  end subroutine check_split_target_dim


  subroutine init_remaining_target_dim(self)
    class(UpscaleHelper), intent(inout) :: self

    integer(i4) :: iCoord
    logical :: is_finalized

    ! check if any target_field is yet uninitialized and try to initialize remaining with help of sourceCoordinate
    do iCoord = 1, size(self%targetCoords)
      ! inquire status of initialization
      is_finalized = self%targetCoords(iCoord)%coord_p%is_finalized()
      if (.not. is_finalized .and. .not. self%targetCoords(iCoord)%coord_p%is_2d()) then
        ! init the remaining target coordinate with help of properties from source coordinate
        call self%targetCoords(iCoord)%coord_p%from_other(self%sourceCoords(iCoord)%coord_p)
      else if (self%targetCoords(iCoord)%coord_p%is_2d() .and. is_finalized .and. &
          self%targetCoords(iCoord)%coord_p%count == defaultCoordCount) then
        call self%targetCoords(iCoord)%coord_p%set_2d_count()
      end if
    end do

  end subroutine init_remaining_target_dim


  subroutine execute_UpscaleHelper(self, data, reshapedMask)
    class(UpscaleHelper), intent(inout) :: self
    real(dp), dimension(:), allocatable, intent(inout) :: data
    logical, dimension(:), allocatable, intent(inout) :: reshapedMask

    !dim (active_cells)
    real(dp), dimension(:), allocatable :: currentArray, currentNanWeights
    logical, dimension(:), allocatable :: currentMask

    integer(i4) :: iCoord, nCoord
    type(CoordUpscaler), pointer :: coord_upscaler => null()

    ! if nothing to be done, simply return data and reshapedMask
    if (.not. self%doUpscale) then
      return
    end if

    ! get the nan-values back in the array, so that the shape of currentArray equals the product of the coordinates
    ! sanity check
    if (product(self%newShapes(:, 1)) /= size(reshapedMask, kind=i8)) then
      log_error(*) "execute_UpscaleHelper: the size of the data (", size(reshapedMask, kind=i8), &
              ") does not match the expected size (", &
      product(self%newShapes(:, 1)), ") - Aborting. Check used Upscaler ", trim(self%name)
      call self%get_stats()
      stop 1
    end if

    if (size(reshapedMask, kind=i8) > int(huge(0_i4), i8)) then
      currentArray = unpack_chunkwise(data, reshapedMask, nodata_dp)
    else
      currentArray = unpack(data, reshapedMask, nodata_dp)
    end if
    call move_alloc(reshapedMask, currentMask)
    deallocate(data)

    ! reorder if not yet in right order
    if (self%doReorder) then
      ! check if the transposing works
      if( product(self%originalShape) /= product(self%newShapes(:, 1))) then
        log_error(*) 'execute_UpscaleHelper: Unexpected configuration, contact developer!'
        call self%get_stats()
        stop 1
      end if
      call reorder_coordinates(currentArray, currentMask, self%originalShape, self%reorderIndices)
    end if

    ! introduce a new array of weights for upscaling multiple coordinates:
    ! it is important to know for each cell that was already upscaled by how many non-nan cells it is supported
    ! this needs to be considered when applying the weights for upscaling
    if (count(self%upscalersIndex /= nodata_i4) > 1_i4) then
      log_debug(*) 'execute_UpscaleHelper: Preparing upscaleHelper "', trim(self%name), &
              '". Activating NanWeights because there is more than one Coordinate upscaler required: ', &
                count(self%upscalersIndex /= nodata_i4)

      allocate(currentNanWeights(size(currentMask, kind=i8)))
      where (currentMask)
        currentNanWeights = 1.0_dp
      else where
        currentNanWeights = 0.0_dp
      end where
    else
      allocate(currentNanWeights(0))
    end if

    ! loop over all coordinates
    nCoord = size(self%newShapes, 1)
    do iCoord = 1, nCoord
      log_trace(*) 'current coordinate: ', self%targetCoords(iCoord)%coord_p%name
      if (.not. self%doCoordUpscale(iCoord)) cycle
      if (self%upscalersIndex(iCoord) /= nodata_i4) then
        ! upscaling is necessary
        ! and link it
        coord_upscaler => MPR_COORD_UPSCALERS(self%upscalersIndex(iCoord))
        log_debug(*) 'Upscaling coordinate: "', trim(coord_upscaler%name), &
                '" at index (', nCoord+1-iCoord, ') using upscale operator: ', &
                trim(self%upscaleOperatorNames(iCoord))

        log_debug(*) "Cells of input array: ", product(self%newShapes(:,iCoord)), &
                " (", self%newShapes(:,iCoord), ")"
        log_debug(*) "Cells of output array: ", product(self%newShapes(:,iCoord+1)), &
                " (", self%newShapes(:,iCoord+1), ")"
        if (size(currentNanWeights, kind=i8) == 0) then
          ! call the upscaling function
          call coord_upscaler%execute(currentArray=currentArray, currentMask=currentMask, &
                  upscaleOperatorName=self%upscaleOperatorNames(iCoord), &
                  currentCoords=self%newShapes(:,iCoord), newCoords=self%newShapes(:,iCoord+1), &
                  indexCoordCurrent=iCoord, indexCoordNew=iCoord)
        else
          call coord_upscaler%execute(currentArray=currentArray, currentMask=currentMask, &
                  currentNanWeights=currentNanWeights, upscaleOperatorName=self%upscaleOperatorNames(iCoord), &
                  currentCoords=self%newShapes(:,iCoord), newCoords=self%newShapes(:,iCoord+1), &
                  indexCoordCurrent=iCoord, indexCoordNew=iCoord)
        end if
      end if
    end do

    ! move the temporary allocation to the intent(out)
    call move_alloc(currentMask, reshapedMask)
    ! pack the temporary data using the updated mask
    data = pack(currentArray, reshapedMask)
    log_subtrace(*) 'execute_UpscaleHelper: upscaled packed data: ', data
    log_subtrace(*) 'execute_UpscaleHelper: upscaled reshapedMask: ', reshapedMask
    if (size(currentNanWeights, kind=i8) /= 0) then
      log_subtrace(*) 'execute_UpscaleHelper: upscaled NanWeights: ', currentNanWeights
    end if

    deallocate(currentArray, currentNanWeights)

    if (self%doSplitTargetCoord) then
      call self%split_target_coordinate()
    end if

  end subroutine execute_UpscaleHelper

  subroutine trim_upscale_operators(self)
    class(UpscaleHelper), intent(inout) :: self
    !< remove item from list of names by given indices
    character(maxNameLength), dimension(:), allocatable :: tempUpscaleOperatorNames

    allocate(tempUpscaleOperatorNames(size(self%upscaleOperatorNames)-1))
    tempUpscaleOperatorNames(1:self%splitTargetIndices(1)) = self%upscaleOperatorNames(1:self%splitTargetIndices(1))
    tempUpscaleOperatorNames(self%splitTargetIndices(1)+1:) = self%upscaleOperatorNames(self%splitTargetIndices(1)+2:)
    call move_alloc(tempUpscaleOperatorNames, self%upscaleOperatorNames)

  end subroutine trim_upscale_operators

  recursive subroutine get_coordUpscaler(sourceCoords, targetCoords, id)
    type(CoordinatePointer), dimension(:), intent(in) :: sourceCoords
    type(CoordinatePointer), dimension(:), intent(in) :: targetCoords
    integer(i4), intent(out) :: id

    type(CoordUpscaler) :: coord_upscaler_
    integer(i4) :: nCoordUpscalers
    character(256) :: coordUpscalerName, coordUpscalerAlias

    ! create the upscaler name
    coordUpscalerName = create_upscaler_name(sourceCoords, targetCoords)
    coordUpscalerAlias = create_upscaler_alias(sourceCoords, targetCoords)

    ! get the number of existing dimupscalers
    nCoordUpscalers = get_n_initialized(MPR_COORD_UPSCALERS)
    ! get the id of the dimupscalers
    id = get_index_in_coord_upscaler(coordUpscalerName, MPR_COORD_UPSCALERS, coordUpscalerAlias)
    if (id > nCoordUpscalers) then
      log_trace(*) 'get_coordUpscaler: Creating CoordUpscaler with name ', trim(coordUpscalerName), &
              ' and alias ', trim(coordUpscalerAlias), ' at id ', id, &
              ' from ', size(sourceCoords), ' source coordinates and ', size(targetCoords), ' target coordinates.'
      ! as this is a recursive subroutine that might init CoordUpscalers in itself
      ! we need to make sure the id is reserved in the global list of CoordUpscalers,
      ! hence we need to init a dummy and non-recursive CoordUpscaler (same name and alias)
      ! and replace that later with the proper one
      coord_upscaler_ = CoordUpscaler(&
              name=coordUpscalerName, &
              id=id, &
              alias=coordUpscalerAlias, &
              initDummyOnly=.true. &
              )
      ! add dummy to the dimupscalers list
      call add_coord_upscaler(coord_upscaler_, id)

      ! create the dimupscaler if it not yet exists
      coord_upscaler_ = CoordUpscaler(&
              name=coordUpscalerName, &
              id=id, &
              alias=coordUpscalerAlias, &
              sourceCoord=sourceCoords, &
              targetCoord=targetCoords &
              )
      ! add it to the dimupscalers list
      log_trace(*) 'get_coordUpscaler: Adding CoordUpscaler with name ', trim(coordUpscalerName), &
              ' and alias ', trim(coordUpscalerAlias), ' at id ', id, &
              ' from ', size(sourceCoords), ' source coordinates and ', size(targetCoords), ' target coordinates.'
      call add_coord_upscaler(coord_upscaler_, id)

    end if

  end subroutine get_coordUpscaler


  subroutine add_upscaler(item, insertIndex)
    class(UpscaleHelper), intent(in) :: item
    integer(i4), intent(in) :: insertIndex

    ! check if there is still room in the vector for one more item
    if (insertIndex >= size(MPR_UPSCALERS)) then
      log_error("(1X,A,A,I0,A,I0)")  "Data array cannot be added at index ('", insertIndex, &
              "'). Vector of data arrays only has length: ", size(MPR_UPSCALERS)
      stop 1
    end if

    ! add item
    MPR_UPSCALERS(insertIndex) = item

  end subroutine add_upscaler

  subroutine split_target_coordinate(self)
    !< splits a target coordinate into its components
    class(UpscaleHelper), intent(inout) :: self

    integer(i4) :: iCoord
    integer(i8) :: indexCounter
    integer(i4), dimension(:), allocatable :: coordIds
    type(CoordinatePointer), dimension(:), allocatable :: tempCoordPointers

    ! two target coordinates were combined for upscaling purposes but now need to be split
    call split_coordinate(self%targetCoords(self%splitTargetIndices(1)), coordIds)
    log_debug(*) 'split_target_coordinate: split target coordinate ', &
            trim(self%targetCoords(self%splitTargetIndices(1))%coord_p%name), ' and got those coordIds: ', coordIds
    if (size(coordIds) /= 2) then
      log_error(*) "The temporary coordinate ", trim(self%targetCoords(self%splitTargetIndices(1))%coord_p%name), &
      " was split back into ", size(coordIds), " coordinates, allowed are only 2."
      stop 1
    end if
    ! allocate the proper target coords
    allocate(tempCoordPointers(size(self%targetCoords) +1))
    ! and set the correct pointers
    indexCounter = 0_i8
    do iCoord=1, size(self%targetCoords)
      if (self%splitTargetIndices(1) == iCoord) then
        ! let the coordinate point to our new coordinate
        call tempCoordPointers(iCoord)%set_coordinate_pointer(coordIds(1))
        call tempCoordPointers(iCoord+1)%set_coordinate_pointer(coordIds(2))
        indexCounter = indexCounter + 1_i8
      else
        tempCoordPointers(iCoord+indexCounter) = self%targetCoords(iCoord)
      end if
    end do
    call move_alloc(tempCoordPointers, self%targetCoords)

  end subroutine split_target_coordinate

  subroutine get_stats_UpscaleHelper(self)
    !< prints some information on the UpscaleHelper
    !< useful for debugging or printing before raising error messages

    !> CoordUpscaler type-bound procedure
    class(UpscaleHelper), intent(in) :: self
    integer(i4) :: i

    log_info("(1X,A,A,A,A,I0)") 'UpscaleHelper "', trim(self%name), '" with id ', self%id
    log_info(*) 'The flag "doUpscale" is: ', self%doUpscale
    if (allocated(self%upscaleOperatorNames)) then
      log_info(*) 'There are ', size(self%upscaleOperatorNames), ' targetCoords: '
      do i=1, size(self%upscaleOperatorNames)
        log_info(*) '  ', i, ': "', trim(self%targetCoords(i)%coord_p%name), '" with sourceCoord "', &
                trim(self%sourceCoords(i)%coord_p%name), '" with upscaleOperatorName "', &
                trim(self%upscaleOperatorNames(i)), '" using CoordUpscaler index', self%upscalersIndex(i)
        log_info(*) '      The flag "doCoordUpscale" is: ', self%doCoordUpscale(i)
      end do
    end if
    log_info(*) 'The flag "doReorder" is: ', self%doReorder
    if (self%doReorder) then
      log_info(*) 'The reorder indices based on original order are: ', self%reorderIndices
      log_info(*) 'The original shape of array (size of coordinates) is: ', self%originalShape
    end if
    log_info(*) 'The flag "doSplitTargetCoord" is: ', self%doSplitTargetCoord
    if (self%doSplitTargetCoord) then
      log_info(*) 'The indices of the (temporary) target coordinate to split into elementary coordinates are: ', &
        self%splitTargetIndices
    end if
    log_info(*) 'The flag "didCombineSourceCoords" is: ', self%didCombineSourceCoords
    if (self%doSplitTargetCoord) then
      log_info(*) 'The indices of (original) sources coordinate combined into 2d coordinate are: ', &
        self%combineSourceIndices
    end if
    if (allocated(self%newShapes)) then
      log_info(*) 'The newShapes array stores the intermediate shapes during upscaling: ', self%newShapes
    end if

  end subroutine get_stats_UpscaleHelper

  subroutine reorder_coordinates(array, mask, shape, order)
    !< reorder the coordinates of an array
    ! TODO: This can be done much faster, I guess: (https://en.wikipedia.org/wiki/In-place_matrix_transposition)
    ! but this is neat and works

    !> the array to be reshaped
    real(dp), dimension(:), intent(inout) :: array
    !> the mask to be reshaped
    logical, dimension(:), intent(inout) :: mask
    !> the shape of the input array (if it would not be packed)
    integer(i8), dimension(:), intent(in) :: shape
    !> the order (indices) of shapes for the reshape (should not be continuous of course)
    !> it can contain values that are bigger than size(order), this refers to new coordinates where broadcasting is
    !>   necessary later on
    integer(i4), dimension(:), intent(in) :: order

    ! helper vectors
    integer(i8), dimension(size(array, kind=i8)) :: arrayIndex, ids, newIndex
    integer(i8), dimension(size(order) + 1) :: reshapedCumProduct
    integer(i8), dimension(size(shape) + 1) :: originalCumProduct
    integer(i4) :: iCoord
    integer(i8) :: i

    ! a vector of ascending ids starting at 0
    ids = [(i-1_i8, i=1_i8, size(array, kind=i8))]
    ! Fortran has 1-based indexing
    arrayIndex = 1_i8

    ! base factor for cumulated product of coordinate sizes
    reshapedCumProduct = 1_i8
    originalCumProduct = 1_i8
    do iCoord = 1, size(shape)
      ! calculate the cumulative product of the size of each coordinate in the original order (source)
      originalCumProduct(iCoord + 1) = originalCumProduct(iCoord) * shape(iCoord)
    end do
    do iCoord = 1, size(order)
      if (order(iCoord) > size(shape)) then
        ! there is a new target coordinate, this is done later by broadcasting, we do not have to care about that now
        reshapedCumProduct(iCoord + 1) = reshapedCumProduct(iCoord)
        cycle
      else
        ! calculate the cumulative product of the size of each coordinate in the future order (target)
        reshapedCumProduct(iCoord + 1) = reshapedCumProduct(iCoord) * shape(order(iCoord))
      end if
      ! ignore size(1) coordinates
      if (shape(order(iCoord)) == 1_i4) then
        cycle
      end if
      ! add the indexing rule for this Coordinate to the indexer
      newIndex = [((mod(ids(i), reshapedCumProduct(iCoord + 1)) / &
              reshapedCumProduct(iCoord)) * originalCumProduct(order(iCoord)), i=1_i8, size(ids, kind=i8))]
      arrayIndex = arrayIndex + newIndex
    end do
    ! reorder by using index vector on current arrays
    array = array(arrayIndex)
    mask = mask(arrayIndex)

  end subroutine reorder_coordinates


  recursive function newCoordUpscaler(name, id, alias, sourceCoord, targetCoord, fromWeightFile, &
            subcellIdsFieldName, weightsFieldName, nSubcellsFieldName, initDummyOnly)
    !< initialization function for type CoordUpscaler
    !< it is recursive as a it can be initialized with multiple targetCoords which first internally initializes
    !< CoordUpscaler for each coordinate pair

    use mo_mpr_global_variables, only: WRITE_WEIGHTS

    character(*), intent(in) :: name
    integer(i4), intent(in) :: id
    character(*), intent(in) :: alias
    type(CoordinatePointer), intent(in), dimension(:), optional :: sourceCoord
    type(CoordinatePointer), intent(in), dimension(:), optional :: targetCoord
    character(*), intent(in), optional :: fromWeightFile
    character(*), intent(in), optional :: subcellIdsFieldName
    character(*), intent(in), optional :: weightsFieldName
    character(*), intent(in), optional :: nSubcellsFieldName
    logical, intent(in), optional :: initDummyOnly
    type(CoordUpscaler) :: newCoordUpscaler

    integer(i4), allocatable, dimension(:) :: ids
    integer(i8), allocatable, dimension(:) :: sizeSourceCoords
    integer(i4) :: iCoord

    ! set the input information
    newCoordUpscaler%name = trim(name)
    newCoordUpscaler%id = id
    newCoordUpscaler%alias = trim(alias)
    newCoordUpscaler%doNeedWeights = .true.
    newCoordUpscaler%mapMethod = defaultMapMethod

    newCoordUpscaler%is_initialized = .true.

    if (present(sourceCoord) .and. present(targetCoord)) then
      ! check for equal lengths
      if (size(sourceCoord) /= size(targetCoord)) then
        log_error("(1X,A,A,A,A,I0,A,I0,A)") 'Error initializing CoordUpscaler "', trim(name), &
                '". Number of source coordinates (', size(sourceCoord), ') and target coordinates (', &
                size(sourceCoord) , ') do not match.'
      end if
      ! create the upscalers for each pair
      if (size(targetCoord, kind=i8) == 1_i8) then
        if (.not. sourceCoord(1_i8)%coord_p%is_2d() .and. .not. targetCoord(1_i8)%coord_p%is_2d()) then
          ! this creates the weights for two 1d coordinate variables
          log_debug (*) 'newCoordUpscaler: init from single pair'
          call newCoordUpscaler%compute_weights_1d(sourceCoord(1_i8)%coord_p, targetCoord(1_i8)%coord_p)
        else if (sourceCoord(1_i8)%coord_p%is_polygon() .and. targetCoord(1_i8)%coord_p%is_polygon()) then
          ! this creates the weights for two polygon-based coordinate variables
          call newCoordUpscaler%compute_weights_poly(sourceCoord(1_i8)%coord_p, targetCoord(1_i8)%coord_p)
        else if (.not. sourceCoord(1_i8)%coord_p%is_polygon() .and. targetCoord(1_i8)%coord_p%is_polygon()) then
          ! this creates the weights for 2d to polygon-based coordinate variables
          call sourceCoord(1_i8)%coord_p%set_polygons_from_2d(centersOnlyArg=.true.)
          call newCoordUpscaler%compute_weights_poly(sourceCoord(1_i8)%coord_p, targetCoord(1_i8)%coord_p)
        else if (sourceCoord(1_i8)%coord_p%is_polygon() .and. .not. targetCoord(1_i8)%coord_p%is_polygon()) then
          ! this creates the weights for polygon-based to 2d coordinate variables
          call targetCoord(1_i8)%coord_p%set_polygons_from_2d()
          call newCoordUpscaler%compute_weights_poly(sourceCoord(1_i8)%coord_p, targetCoord(1_i8)%coord_p)
        else if (.not. sourceCoord(1_i8)%coord_p%is_polygon() .and. .not. targetCoord(1_i8)%coord_p%is_polygon()) then
          ! this creates the weights for two 2d coordinate variables
          call newCoordUpscaler%compute_weights_2d(sourceCoord(1_i8)%coord_p, targetCoord(1_i8)%coord_p)
        else
          log_error(*) 'unexpected error, contact developer'
          stop 1
        end if
      else
        log_debug (*) 'newCoordUpscaler: init from ', size(targetCoord), 'pairs'
        allocate(ids(size(targetCoord)), sizeSourceCoords(size(sourceCoord)))
        ! create upscaler for each coordinate pair
        do iCoord=1, size(targetCoord)
          call get_coordUpscaler([sourceCoord(iCoord)], [targetCoord(iCoord)], ids(iCoord))
        end do
        ! log the size of the source coordinate, this is needed for getting the indices right during combination
        do iCoord=1, size(sourceCoord)
          sizeSourceCoords(iCoord) = sourceCoord(iCoord)%coord_p%count
        end do
        call newCoordUpscaler%combine_weights(ids, sizeSourceCoords)
      end if

    else if (present(fromWeightFile)) then
      ! this reads the weights
      call newCoordUpscaler%read_weights(fromWeightFile, subcellIdsFieldName, weightsFieldName, nSubcellsFieldName)
    else if (.not. present(initDummyOnly)) then
      log_error(*) 'Error initializing CoordUpscaler "', trim(name), '". Either pass arguments',  &
              '"sourceCoord" and "targetCoord" or "fromWeightFile", "subcellIdsFieldName",', &
              ' "weightsFieldName" and "nSubcellsFieldName" or "initDummyOnly".'
      stop 1
    end if

    if (.not. present(initDummyOnly) .and. WRITE_WEIGHTS) call newCoordUpscaler%write_weights()

  end function newCoordUpscaler

  function is_finalized_CoordUpscaler(self) result(is_finalized)
    class(CoordUpscaler), intent(in) :: self
    logical :: is_finalized

    is_finalized = allocated(self%weights)
  end function is_finalized_CoordUpscaler

  subroutine reset_CoordUpscaler(self)
    class(CoordUpscaler), intent(inout) :: self

    if (self%is_initialized) then
      self%name = ''
      self%id=0
      self%is_initialized = .false.
      self%alias = ''
      if (allocated(self%subcells)) deallocate(self%subcells)
      if (allocated(self%ids)) deallocate(self%ids)
      if (allocated(self%weights)) deallocate(self%weights)
    end if

  end subroutine reset_CoordUpscaler

  subroutine execute_CoordUpscaler(self, currentArray, currentMask, currentNanWeights, upscaleOperatorName, &
                currentCoords, newCoords, indexCoordCurrent, indexCoordNew)
    !< this subroutine requires three packed input arrays
    !< with the values, their mask and the weights for each value

    !> the procedure is attributed to the class CoordUpscaler
    class(CoordUpscaler), intent(inout) :: self
    !dim of the following dummy variables: packed cells of source array
    !> packed (no mask used) values of source array
    real(dp), dimension(:), allocatable, intent(inout) :: currentArray
    !> mask declaring whether all subcells have been nan
    logical, dimension(:), allocatable, intent(inout) :: currentMask
    !> weights for each cell determining how many non-nan cells were used to compute the current cell
    !> 0: none 1: only nan
    real(dp), dimension(:), allocatable, intent(inout), optional :: currentNanWeights
    !> name of upscaling operator to be used for referencing function
    character(*), intent(in) :: upscaleOperatorName
    !> vector of all dimension sizes for the input array
    integer(i8), dimension(:), intent(in) :: currentCoords
    !> vector of all dimension sizes for the output array
    integer(i8), dimension(:), intent(in) :: newCoords
    !> F-style index of the first coordinate to include in upscaling
    integer(i4), intent(in) :: indexCoordCurrent
    !> F-style index of the last coordinate to include in upscaling
    integer(i4), intent(in) :: indexCoordNew

    !dim of the following dummy variables: packed cells of target array
    real(dp), dimension(:), allocatable :: newArray
    logical, dimension(:), allocatable :: newMask
    real(dp), dimension(:), allocatable :: newNanWeights

    integer(i8), dimension(:), allocatable :: effectiveSubcellIds, subcellIds, idAddon
    procedure(upscale_func_alias), pointer :: upscaleFunc => null()
    real(dp) :: p
    integer(i8) :: iCell, idMultiplier
    integer(i4) :: nCoord

    ! this check is done for case, when weights are read from external source,
    ! applied to coordinates by name only, without checking of coordinate bounds of coordinate actually existing
    nCoord = size(newCoords)
    log_trace(*) 'execute_CoordUpscaler: received currentCoords:', currentCoords
    log_trace(*) 'execute_CoordUpscaler: received newCoords:', newCoords
    log_trace(*) 'execute_CoordUpscaler: received indexCoordCurrent:', indexCoordCurrent
    log_trace(*) 'execute_CoordUpscaler: received indexCoordNew:', indexCoordNew
    if (product(newCoords(nCoord+1_i4-indexCoordNew:nCoord+1_i4-indexCoordCurrent)) /= size(self%subcells, kind=i8)) then
      log_error('(1X,A,A,A,A,I0,A,I0,A,A)')  "execute_CoordUpscaler: Coordinate upscaler '", trim(self%alias), &
              "' cannot be executed as the size of the n_subcells array (", &
              size(self%subcells, kind=i8), &
              ") and the size of the target array (slice) (", newCoords(indexCoordNew), ") do not match.", &
              " Check, if the weights file is correct, if initialized by file."
      stop 1
    end if

    allocate(subcellIds(product(newCoords)))
    allocate(idAddon(product(newCoords)))

    ! this computes the helper variables (subcellIds, idMultiplier, idAddon)
    ! used to compute the effectiveSubcellIds
    call compute_index_helpers(currentCoords, newCoords, indexCoordCurrent, indexCoordNew, idAddon, subcellIds, idMultiplier)
    ! this translates a string for the upscaler into a pointer to procedure and optionally also yielding
    ! p-norm value
    call get_upscale_func(upscaleOperatorName, upscaleFunc, p, self%doNeedWeights .or. present(currentNanWeights))

    ! prepare the new arrays
    allocate(newArray(product(newCoords)))
    allocate(newMask(product(newCoords)))

    if (present(currentNanWeights)) then
      allocate(newNanWeights(product(newCoords)))
      ! loop over cells
      do iCell=1, size(newArray, kind=i8)
        if (self%subcells(subcellIds(iCell)) == 0_i8) then
          newArray(iCell) = nodata_dp
          newMask(iCell) = .false.
          newNanWeights(iCell) = 0.0_dp
        else
          ! the effective index array for the packed currentArray declaring
          ! which subcells are needed for iCell of newArray
          effectiveSubcellIds = (self%ids(1:self%subcells(subcellIds(iCell)), &
                                          subcellIds(iCell)) - 1) * &
                  idMultiplier + 1 + idAddon(iCell)
          log_subtrace(*) 'execute_CoordUpscaler: at newCell (', iCell, ') using the ids (', &
                  effectiveSubcellIds, ') with these values (', currentArray(effectiveSubcellIds), &
                  ') and those weights (', self%weights(1:self%subcells(subcellIds(iCell)), subcellIds(iCell)), &
                  ') in mode with doNeedWeights and with NanWeights'
          ! now execute the upscaling for the current cell
          if (mod(iCell, maxval([size(newArray, kind=i8) / 10_i8, 1_i8])) == 0_i8 .and. &
                  size(newArray, kind=i8) > 1000000_i8) then
            log_debug(*) 'execute_CoordUpscaler: Upscaling iCell: ', iCell, '/', size(newArray, kind=i8)
          end if
          ! now execute the upscaling for the current cell
          call wrap_weighted_upscale(&
                  upscaleFunc, &
                  p, &
                  self%weights(1:self%subcells(subcellIds(iCell)), subcellIds(iCell)), &
                  ! intent in
                  currentArray(effectiveSubcellIds), &
                  currentMask(effectiveSubcellIds), &
                  currentNanWeights(effectiveSubcellIds), &
                  ! intent out
                  newArray(iCell), &
                  newMask(iCell), &
                  newNanWeights(iCell) &
                  )
        end if
      end do
      call move_alloc(newNanWeights, currentNanWeights)
    else
      if (self%doNeedWeights) then
        ! loop over cells
        do iCell=1, size(newArray, kind=i8)
          if (self%subcells(subcellIds(iCell)) == 0_i8) then
            newArray(iCell) = nodata_dp
            newMask(iCell) = .false.
          else
            ! the effective index array for the packed currentArray declaring
            ! which subcells are needed for iCell of newArray
            effectiveSubcellIds = (self%ids(1:self%subcells(subcellIds(iCell)), subcellIds(iCell)) - 1) * &
                    idMultiplier + 1 + idAddon(iCell)
            log_subtrace(*) 'execute_CoordUpscaler: at newCell (', iCell, ') using the ids (', &
                    effectiveSubcellIds, ') with these values (', currentArray(effectiveSubcellIds), &
                    ') and those weights (', self%weights(1:self%subcells(subcellIds(iCell)), subcellIds(iCell)), &
                    ') in mode with doNeedWeights and without NanWeights'
            ! now execute the upscaling for the current cell
            if (mod(iCell, maxval([size(newArray, kind=i8) / 10_i8, 1_i8])) == 0_i8 .and. &
                    size(newArray, kind=i8) > 1000000_i8) then
              log_debug(*) 'execute_CoordUpscaler: Upscaling iCell: ', iCell, '/', size(newArray, kind=i8)
            end if
            call wrap_weighted_upscale(&
                    func=upscaleFunc, &
                    p=p, &
                    weightsIn=self%weights(1:self%subcells(subcellIds(iCell)), subcellIds(iCell)), &
                    ! intent in
                    sliceIn=currentArray(effectiveSubcellIds), &
                    maskIn=currentMask(effectiveSubcellIds), &
                    ! intent out
                    valueOut=newArray(iCell), &
                    maskOut=newMask(iCell) &
                    )
          end if

        end do
      else
        ! loop over cells
        do iCell=1, size(newArray, kind=i8)
          if (self%subcells(subcellIds(iCell)) == 0_i8) then
            newArray(iCell) = nodata_dp
            newMask(iCell) = .false.
          else
            ! the effective index array for the packed currentArray declaring which subcells are needed for iCell of newArray
            effectiveSubcellIds = (self%ids(1:self%subcells(subcellIds(iCell)), subcellIds(iCell)) - 1) * &
                    idMultiplier + 1 + idAddon(iCell)
            log_subtrace(*) 'execute_CoordUpscaler: at newCell (', iCell, ') using the ids (', &
                    effectiveSubcellIds, ') with these values (', currentArray(effectiveSubcellIds), &
                    ') and those weights (', self%weights(1:self%subcells(subcellIds(iCell)), subcellIds(iCell)), &
                    ') in mode without doNeedWeights and without NanWeights'
            ! now execute the upscaling for the current cell
            if (mod(iCell, maxval([size(newArray, kind=i8) / 10_i8, 1_i8])) == 0_i8 .and. &
                size(newArray, kind=i8) > 1000000_i8) then
              log_debug(*) 'execute_CoordUpscaler: Upscaling iCell: ', iCell, '/', size(newArray, kind=i8)
            end if
            ! now execute the upscaling for the current cell
            call wrap_upscale(&
                    upscaleFunc, &
                    p, &
                    ! intent in
                    sliceIn=currentArray(effectiveSubcellIds), &
                    maskIn=currentMask(effectiveSubcellIds), &
                    ! intent out
                    valueOut=newArray(iCell), &
                    maskOut=newMask(iCell) &
                    )
          end if
        end do
      end if
    end if

    ! move the allocation from the new to the current array (implicitly deallocating the new...)
    call move_alloc(newArray, currentArray)
    call move_alloc(newMask, currentMask)
    deallocate(subcellIds, idAddon)

  end subroutine execute_CoordUpscaler

  subroutine compute_index_helpers(currentCoords, newCoords, indexCoordCurrent, indexCoordNew, &
          idAddon, subcellIds, idMultiplier)
    !< computes the index helpers for upscaling a coordinate
    !< the index helpers are required to select the correct ids for upscaling from the packed input array
    integer(i8), dimension(:), intent(in) :: currentCoords
    integer(i8), dimension(:), intent(in) :: newCoords
    integer(i4), intent(in) :: indexCoordCurrent
    integer(i4), intent(in) :: indexCoordNew
    integer(i8), dimension(:), intent(inout) :: idAddon
    integer(i8), dimension(:), intent(inout), optional :: subcellIds
    integer(i8), intent(out), optional :: idMultiplier

    integer(i8), dimension(size(idAddon, kind=i8)) :: cellIds
    integer(i8) :: iCell, nCoord, a, b

    ! get the size of the coordinates
    nCoord = size(currentCoords, kind=i8)
    ! get an ascending index for each packed target cell
    cellIds = [(iCell, iCell=1, size(idAddon, kind=i8))]
    ! this is the multiplier, that is later on used as a factor for the subcellIds directly
    a = product(currentCoords(nCoord+1-indexCoordCurrent+1:))
    ! generating the value for selecting the correct subcellIds index (last dim in self%subcellIds)
    b = product(newCoords(nCoord+1-indexCoordNew: nCoord+1-indexCoordCurrent))
    log_trace(*) 'compute_index_helpers: a,b for subcellIds', a, b

    if (present(idMultiplier)) then
      idMultiplier = a
      log_trace(*) 'compute_index_helpers: idMultiplier', idMultiplier
    end if

    if (present(subcellIds)) then
      subcellIds = permutate(cellIds, a, b) + 1
      log_subtrace(*) 'compute_index_helpers: subcellIds', subcellIds
    end if

    ! this is the addon, that is later on used on the subcellIds
    idAddon = 0_i8
    if (a > 1_i8) then
      idAddon = idAddon + permutate(cellIds, 1_i8, a) * 1_i8
    end if
    if (a*b < product(newCoords(:))) then
      idAddon = idAddon + permutate(cellIds, a*b, product(newCoords(:)) / (a*b)) * &
              product(currentCoords(nCoord+1-indexCoordNew:))
    end if
    log_subtrace(*) 'compute_index_helpers: idAddon', idAddon

  end subroutine compute_index_helpers

  elemental function permutate(ids, a, b)
    !< little helper function creating a permutation of indices
    integer(i8), intent(in) :: ids
    integer(i8), intent(in) :: a
    integer(i8), intent(in) :: b
    integer(i8) :: permutate

    permutate = mod(ids-1, a*b)/a

  end function permutate

  subroutine add_coord_upscaler(item, insertIndex)
    class(CoordUpscaler), intent(in) :: item
    integer(i4), intent(in) :: insertIndex

    ! check if there is still room in the vector for one more item
    if (insertIndex >= size(MPR_COORD_UPSCALERS)) then
      log_error("(1X,A,A,A,A,I0,A,I0)")  "Coordinate upscaler '", trim(item%name), "'cannot be added at index ('", &
              insertIndex, "'). Vector of coordinate upscalers only has length: ", size(MPR_COORD_UPSCALERS)
      stop 1
    end if

    ! add item
    MPR_COORD_UPSCALERS(insertIndex) = item

  end subroutine add_coord_upscaler

  subroutine compute_weights_1d(self, sourceCoord, targetCoord)
    !< compute the weights for two regular 1d coordinate variables
    class(CoordUpscaler), intent(inout) :: self
    type(Coordinate), intent(in) :: sourceCoord
    type(Coordinate), intent(in) :: targetCoord
    ! pointer to the function to use to get lower and upper source coordinate index of current target cell
    procedure(get_source_index_alias), pointer :: get_lower_source_index_func => null()
    procedure(get_source_index_alias), pointer :: get_upper_source_index_func => null()
    integer(i8) :: cachedIndex
    integer(i4) :: caseId
    real(dp) :: targetStep, sourceStep

    self%mapMethod = 'Conservative remapping'

    ! check if other_dim fully covers targetCoord
    call check_within_bounds(targetCoord, sourceCoord%bounds)

    log_debug(*) 'compute_weights_1d: Computing weights for CoordUpscaler: ', self%name
    log_trace(*) "compute_weights_1d: targetCoord%bounds", targetCoord%bounds
    log_trace(*) "compute_weights_1d: sourceCoord%bounds", sourceCoord%bounds

    ! handle the different cases for constant step sizes and ascending order of values
    caseId = 0_i4
    if (ne(targetCoord%step, nodata_dp)) then
      ! constant target coord step size
      caseId = caseId + 1_i4
      targetStep = targetCoord%step
    end if
    if (ne(sourceCoord%step, nodata_dp)) then
      ! constant source coord step size
      caseId = caseId + 2_i4
      sourceStep = sourceCoord%step
      if (sourceCoord%is_ascending()) then
        ! ascending order of source coord values
        log_trace(*) "compute_weights_1d: using get_source_index_asc_eq"
        get_upper_source_index_func => get_upper_source_index_asc_eq
        get_lower_source_index_func => get_lower_source_index_asc_eq
        cachedIndex = 1_i8
      else
        ! descending order of source coord values
        log_trace(*) "compute_weights_1d: using get_source_index_desc_eq"
        get_upper_source_index_func => get_upper_source_index_desc_eq
        get_lower_source_index_func => get_lower_source_index_desc_eq
        cachedIndex = sourceCoord%count
      end if
    else
      if (sourceCoord%is_ascending()) then
        ! ascending order of source coord values (source_un_asc)
        log_trace(*) "compute_weights_1d: using get_source_index_asc_un"
        get_upper_source_index_func => get_upper_source_index_asc_un
        get_lower_source_index_func => get_lower_source_index_asc_un
        cachedIndex = 1_i8
      else
        ! descending order of source coord values (source_un_desc)
        log_trace(*) "compute_weights_1d: using get_source_index_desc_un"
        get_upper_source_index_func => get_upper_source_index_desc_un
        get_lower_source_index_func => get_lower_source_index_desc_un
        cachedIndex = sourceCoord%count
      end if
    end if

    select case(caseId)
    case(3)
      ! check if the constant step sizes are multiples, then init the values in one go (former mHM-MPR-style)
      call self%check_step_multiples(sourceCoord, targetCoord)
      if (self%doNeedWeights) then
        ! initialize our target arrays (there are no multiples)
        ! get the ratio of both, this is also the max no of self%subcells
        call self%init_values(int(ceiling(targetStep / sourceStep), kind=i8), targetCoord%count)
        call self%compute_weights_1d_(sourceCoord, targetCoord, &
                get_lower_source_index_func, get_upper_source_index_func, cachedIndex)
      end if
    case(0 : 2)
      ! initialize our target arrays, guess the subcell dimension size
      call self%init_values(int(sourceCoord%count/targetCoord%count*1.5_dp, kind=i8), targetCoord%count)
      call self%compute_weights_1d_(sourceCoord, targetCoord, &
                get_lower_source_index_func, get_upper_source_index_func, cachedIndex)
      ! adapt the sizes
      call resize_weights_array(self%weights, self%ids, self%subcells)
    end select
    log_debug(*) 'compute_weights_1d: first ', minval([targetCoord%count, 10_i8]), ' subcells: ', &
            self%subcells(1:minval([targetCoord%count, 10_i8]))
    log_debug(*) 'compute_weights_1d: first ', minval([self%subcells(1), 10_i8]), ' ids of first target cell: ', &
            self%ids(1:minval([self%subcells(1), 10_i8]), 1_i8)

  end subroutine compute_weights_1d

  subroutine check_step_multiples(self, sourceCoord, targetCoord)
    !< check if source and target Coordinates' step sizes are multiples, then set values accordingly
    !< assume, two 1d Coordinates with a constant step size
    class(CoordUpscaler), intent(inout) :: self
    type(Coordinate), intent(in) :: sourceCoord
    type(Coordinate), intent(in) :: targetCoord

    real(dp) :: ratio1, ratio2, offset
    logical :: conditionDownscale, conditionUpscale, conditionCommon

    self%doNeedWeights = .true.
    ! set some ratios of the step widths
    ratio1 = abs(sourceCoord%step/targetCoord%step)
    ratio2 = abs(targetCoord%step/sourceCoord%step)
    ! offset of source lower bound from target lower bound
    offset = abs(minval(sourceCoord%bounds) - minval(targetCoord%bounds)) / maxval([ratio1, ratio2])

    ! check if ratio and offset rate equals an integer value and offset is a multiple of either
    ! TODO: how to reformulate compute_weights_1d_multiple, so it also works for cases where only first
    !   part of condition is met (test case 3 in test_compute_weights_1d)
    conditionCommon = (abs(anint(offset) - offset) < maxTolerance) .and. &
            ((abs(anint(offset/sourceCoord%step) - offset/sourceCoord%step) < maxTolerance) .or. &
             (abs(anint(offset/targetCoord%step) - offset/targetCoord%step) < maxTolerance))
    conditionDownscale = eq(anint(ratio1), ratio1) .and. conditionCommon
    conditionUpscale = eq(anint(ratio2), ratio2) .and. conditionCommon

    if (conditionDownscale) then
      ! downscaling, one subcell only for each target cell
      log_debug(*) "check_step_multiples: found source Coordinate '", trim(sourceCoord%name), &
              "' (step: ", sourceCoord%step, ") to be a multiple (", ratio1, &
              ") of the target Coordinate '", trim(targetCoord%name), "' (step: ", targetCoord%step, &
              ") with aligned offset: ", offset
      call self%compute_weights_1d_multiple(ratio1, ratio2, int(offset, kind=i8), sourceCoord, targetCoord)
    else if (conditionUpscale) then
      ! traditional upscaling mHM-MPR-style (L1 is multiple of L0 resolution)
      log_debug(*) "check_step_multiples: found target Coordinate '", trim(targetCoord%name), &
        "' (step: ", targetCoord%step, ") to be a multiple (", ratio2, &
        ") of the source Coordinate '", trim(sourceCoord%name), "' (step: ", sourceCoord%step, &
        ") with aligned offset: ", offset
      call self%compute_weights_1d_multiple(ratio1, ratio2, int(offset, kind=i8), sourceCoord, targetCoord)
    end if

  end subroutine check_step_multiples

  subroutine compute_weights_1d_multiple(self, ratio1, ratio2, offset, sourceCoord, targetCoord)
    class(CoordUpscaler), intent(inout) :: self
    real(dp), intent(in) :: ratio1
    real(dp), intent(in) :: ratio2
    integer(i8), intent(in) :: offset
    type(Coordinate), intent(in) :: sourceCoord
    type(Coordinate), intent(in) :: targetCoord

    integer(i8) :: a, b, high, j
    integer(i8), dimension(:), allocatable :: ids

    ! coefficients used later on for id creation
    a = maxval([int(ratio1, kind=i8), 1_i8])
    b = maxval([int(ratio2, kind=i8), 1_i8])
    high = targetCoord%count * b
    log_trace(*) 'ratio1 ', ratio1, ', ratio2 ', ratio2, ', offset ', offset, ', a ', a, ', b ', b, ', high ', high

    if (sourceCoord%is_ascending() .eqv. targetCoord%is_ascending()) then
      ids = permutate([(j, j=1_i8, high, 1_i8)], a, high/a) + 1_i8 + offset
    else
      ids = permutate([(j, j=high, 1_i8, -1_i8)], a, high/a) + 1_i8 + offset
    end if

    allocate(self%subcells(targetCoord%count))
    allocate(self%ids(b, targetCoord%count))
    allocate(self%weights(b, targetCoord%count))

    self%subcells(:) = b
    self%weights(:, :) = minval([ratio1, 1.0_dp])
    self%ids(:, :) = reshape(ids, [b, high/b])
    self%doNeedWeights = .false.

    deallocate(ids)


  end subroutine compute_weights_1d_multiple

  subroutine compute_weights_2d(self, sourceCoord, targetCoord)
    !< compute the weights for two regular 2d coordinate variables by
    !< creating a CoordUpscaler for each 1d coordinate pair, combining them and making a copy
    !< of the ids, subcells and weights for the current CoordUpscaler
    class(CoordUpscaler), intent(inout) :: self
    type(Coordinate), intent(in) :: sourceCoord
    type(Coordinate), intent(in) :: targetCoord

    type(CoordinatePointer), dimension(:), allocatable :: sourceCoords
    type(CoordinatePointer), dimension(:), allocatable :: targetCoords
    integer(i4) :: id, iCoord

    ! there should already have been checks on whether the subDims are already correctly paired
    ! init the pointers

    allocate(sourceCoords(size(sourceCoord%subDims)), targetCoords(size(sourceCoord%subDims)))
    do iCoord=1, size(sourceCoords)
      call sourceCoords(iCoord)%set_coordinate_pointer(sourceCoord%subDims(iCoord))
    end do
    do iCoord=1, size(targetCoords)
      call targetCoords(iCoord)%set_coordinate_pointer(targetCoord%subDims(iCoord))
    end do
    ! create the upscalers based on the the subDims
    call get_coordUpscaler(sourceCoords, targetCoords, id)
    ! use the properties of the just created UpscaleHelper
    log_debug(*) 'creating upscaler for two 2d-coordinate variables using upscaler ', trim(mpr_coord_upscalers(id)%name)
    call self%from_other(mpr_coord_upscalers(id))
    ! adapt the sizes
    call resize_weights_array(self%weights, self%ids, self%subcells)

  end subroutine compute_weights_2d

  subroutine compute_weights_poly(self, sourceCoord, targetCoord)
    !< compute the weights for remapping of two polygon coordinates (unstructured grids)
    !< the algorithm currently uses a simple ray casting algorithm
    class(CoordUpscaler), intent(inout) :: self
    type(Coordinate), intent(in) :: sourceCoord
    type(Coordinate), intent(in) :: targetCoord

    integer(i8) :: iSubcell, iPolygon, currentPolygon, lastiPolygon, maxSubcells, nNodes
    integer(i4) :: result
    real(dp), dimension(2) :: point
    real(dp), dimension(4) :: boundingBox
    real(dp), dimension(:,:), allocatable :: polygonNodes
    integer(i4) :: orientation
    integer(i8), dimension(:), allocatable :: indices, insertIndices
    logical :: useBoundingBox

    log_debug(*) 'compute_weights_poly: computing weights for 2d-coordinate to polygon coordinate variables ', &
            trim(sourceCoord%name), ' to ', trim(targetCoord%name)
    self%mapMethod = 'Nearest destination to source'

    useBoundingBox = .true.
    allocate(polygonNodes(maxval(targetCoord%nodes), 2))

    do iPolygon=1, targetCoord%count
      nNodes = targetCoord%nodes(iPolygon)
      ! select correct polygonNodes
      polygonNodes(1:nNodes,1) = targetCoord%cornersCoord1(1:nNodes, iPolygon)
      polygonNodes(1:nNodes,2) = targetCoord%cornersCoord2(1:nNodes, iPolygon)
      orientation = orientpoly(polygonNodes(1:nNodes,:))
      if (orientation /= -1_i4) then
        useBoundingBox = .false.
        exit
      end if
    end do

    ! restrain points to be considered for search by min and max for x and y
    if (useBoundingBox) then
      if (all(targetCoord%nodes == targetCoord%nodes(1))) then
        ! all polygons have same number of corners
        boundingBox(1) = minval(targetCoord%cornersCoord1(1:targetCoord%nodes(1), :)) - eps_dp
        boundingBox(2) = maxval(targetCoord%cornersCoord1(1:targetCoord%nodes(1), :)) + eps_dp
        boundingBox(3) = minval(targetCoord%cornersCoord2(1:targetCoord%nodes(1), :)) - eps_dp
        boundingBox(4) = maxval(targetCoord%cornersCoord2(1:targetCoord%nodes(1), :)) + eps_dp
        log_trace(*) 'compute_weights_poly: all polygons have same number of nodes: ', targetCoord%nodes(1)
       else
        ! polygons have different numbers of corners
        boundingBox = [targetCoord%cornersCoord1(1, 1), targetCoord%cornersCoord1(1, 1), &
                       targetCoord%cornersCoord2(1, 1), targetCoord%cornersCoord2(1, 1)]
        do iPolygon=1, targetCoord%count
          boundingBox(1) = minval([boundingBox(1), minval(targetCoord%cornersCoord1(1:targetCoord%nodes(iPolygon), iPolygon)) - eps_dp])
          boundingBox(2) = maxval([boundingBox(2), maxval(targetCoord%cornersCoord1(1:targetCoord%nodes(iPolygon), iPolygon)) + eps_dp])
          boundingBox(3) = minval([boundingBox(3), minval(targetCoord%cornersCoord2(1:targetCoord%nodes(iPolygon), iPolygon)) - eps_dp])
          boundingBox(4) = maxval([boundingBox(4), maxval(targetCoord%cornersCoord2(1:targetCoord%nodes(iPolygon), iPolygon)) + eps_dp])
        end do
        log_trace(*) 'compute_weights_poly: polygons have different number of nodes: ', maxval(targetCoord%nodes)
      end if
    end if

    allocate(indices(sourceCoord%count))
    indices = 0_i8
    lastiPolygon = 1_i8
    orientation = -1_i4
    allocate(self%subcells(targetCoord%count))
    self%subcells = 0_i8

    ! loop over source cells
    do iSubcell=1, sourceCoord%count
      ! check if within bounds
      if (useBoundingBox) then
        if ((sourceCoord%centersCoord1(iSubcell) < boundingBox(1)) &
                .or. (sourceCoord%centersCoord1(iSubcell) > boundingBox(2)) &
                .or. (sourceCoord%centersCoord2(iSubcell) < boundingBox(3)) &
                .or. (sourceCoord%centersCoord2(iSubcell) > boundingBox(4))) then
          log_trace(*) 'compute_weights_poly: point ', &
                  [sourceCoord%centersCoord1(iSubcell), sourceCoord%centersCoord2(iSubcell)], &
          ' not in boundingBox: ', boundingBox
          cycle
        end if
      end if

      do iPolygon=1, targetCoord%count
        point = [sourceCoord%centersCoord1(iSubcell), sourceCoord%centersCoord2(iSubcell)]
        ! check the index of last polygon first
        currentPolygon = mod(lastiPolygon + iPolygon - 2, targetCoord%count) + 1_i8
        ! do not recalculate if same polygon is used
        if (iSubcell == 1_i8 .or. currentPolygon /= lastiPolygon) then
          nNodes = targetCoord%nodes(currentPolygon)
          ! select correct polygonNodes
          polygonNodes(1:nNodes,1) = targetCoord%cornersCoord1(1:nNodes, currentPolygon)
          polygonNodes(1:nNodes,2) = targetCoord%cornersCoord2(1:nNodes, currentPolygon)
          orientation = orientpoly(polygonNodes(1:nNodes,:))
          log_trace(*) 'compute_weights_poly: got orientation ', orientation, ' for polygon ', currentPolygon, &
                  ' with nodes ', polygonNodes(1:nNodes,:), ' for point ', point

          if (orientation == 0_i4) then
            ! the polygon covers a pole, so we need to modify the coordinates
            ! warning: this only applies to polygons whose points all have the same latitude
            polygonNodes(1:nNodes,:) = mod_pole(polygonNodes(1:nNodes,:))
          else if (orientation == 1_i4) then
            ! the polygon is clockwise (this happens, if polygon crosses 180 meridian), shift coordinates
            polygonNodes(1:nNodes,1) = mod_shift(polygonNodes(1:nNodes,1))
          end if
        end if
        if (orientation == 1_i4) then
          point(1) = mod_shift(point(1))
        end if
        call inpoly(point, polygonNodes(1:nNodes,:), result)
        if (result > 0_i4) then
          log_trace(*) 'compute_weights_poly: subcell ', iSubcell, ' in target polygon ', currentPolygon
          indices(iSubcell) = currentPolygon
          ! advance counter
          self%subcells(currentPolygon) = self%subcells(currentPolygon) + 1_i8
          exit
        end if
      end do
      lastiPolygon = currentPolygon
    end do
    maxSubcells = maxval(self%subcells)
    log_debug(*) 'compute_weights_poly: setting ids and weights for Upscaler with dimensions (', &
        maxSubcells, 'x', targetCoord%count, ')'

    ! set the ids
    allocate(self%ids(maxSubcells, targetCoord%count))
    self%ids = 0_i8
    allocate(insertIndices(targetCoord%count))
    insertIndices = 1_i8
    do iSubcell=1, sourceCoord%count
      if (indices(iSubcell) > 0_i8) then
        self%ids(insertIndices(indices(iSubcell)), indices(iSubcell)) = iSubcell
        insertIndices(indices(iSubcell)) = insertIndices(indices(iSubcell)) + 1_i8
      end if
    end do
    deallocate(insertIndices, indices, polygonNodes)

    ! set the weights
    allocate(self%weights(maxSubcells, targetCoord%count))
    self%weights = 0.0_dp
    do iPolygon=1, targetCoord%count
      if (self%subcells(iPolygon) > 0_i8) then
        log_trace(*) 'compute_weights_poly: Found ', self%subcells(iPolygon), ' subcells in target polygon ', iPolygon
        self%weights(1:self%subcells(iPolygon), iPolygon) = 1.0_dp / self%subcells(iPolygon)
      end if
    end do

  end subroutine compute_weights_poly

  subroutine compute_weights_1d_(self, sourceCoord, targetCoord, &
          get_lower_source_index_func, get_upper_source_index_func, cachedIndex)
    !< this routine generates the weights for the case where both coords have an equal step size over
    !< their respective cells
    class(CoordUpscaler), intent(inout) :: self
    type(Coordinate), intent(in) :: sourceCoord
    type(Coordinate), intent(in) :: targetCoord
    procedure(get_source_index_alias), intent(in), pointer :: get_lower_source_index_func
    procedure(get_source_index_alias), intent(in), pointer :: get_upper_source_index_func
    integer(i8), intent(in) :: cachedIndex

    integer(i8) :: i, j, lowerSourceIndex, upperSourceIndex, increment, ascFlag, offset1, offset2
    real(dp) :: lowerTargetBound, upperTargetBound

    ! init the upperSourceIndex
    upperSourceIndex = cachedIndex
    increment = 1_i8
    offset1 = -1_i8
    offset2 = 0_i8
    if (.not. sourceCoord%is_ascending()) then
      increment = -1_i8
      offset1 = 0_i8
      offset2 = -1_i8
    end if
    ascFlag = 0_i8
    if (.not. targetCoord%is_ascending()) then
      ascFlag = 1_i8
    end if
    ! loop over targetCoord cell (ends), they are 0-based
    do i = 1, targetCoord%count
      ! get the lower and upper sourceIndex:
      ! compare cell's lower and upper bound (lowerTargetBound and upperTargetBound)
      ! with sourceCoords from currentSourceIndex onwards using step
      lowerTargetBound = targetCoord%values(i-1+ascFlag)
      upperTargetBound = targetCoord%values(i-ascFlag)
      log_trace(*) 'compute_weights_1d_: getting lowerSourceIndex with targetBound:', &
              lowerTargetBound, ', values: ', sourceCoord%values(:), ' (step: ', sourceCoord%step, &
              ') and initial index: ', upperSourceIndex
      lowerSourceIndex = get_lower_source_index_func(lowerTargetBound, sourceCoord%values(:), &
               upperSourceIndex, sourceCoord%step)
      log_trace(*) 'compute_weights_1d_: getting upperSourceIndex with targetBound:', &
              upperTargetBound, ', values: ', sourceCoord%values(:), ' (step: ', sourceCoord%step, &
              ') and initial index: ', lowerSourceIndex
      upperSourceIndex = get_upper_source_index_func(upperTargetBound, sourceCoord%values(:), &
              lowerSourceIndex, sourceCoord%step)
      ! get subcells
      self%subcells(i) = abs(upperSourceIndex - lowerSourceIndex) + 1
      log_trace(*) 'compute_weights_1d_: for cell: ', i, &
              ' got lowerSourceIndex: ', lowerSourceIndex, ' and upperSourceIndex:', upperSourceIndex, &
              ' with increment:', increment
      ! recalculate target arrays, if needed
      if (self%subcells(i) > size(self%ids, 1, kind=i8)) then
        call resize_weights_array(self%weights, self%ids, self%subcells, &
                self%subcells(i), targetCoord%count)
      end if
      ! get ids
      self%ids(1: self%subcells(i), i) = [(j, j = lowerSourceIndex, upperSourceIndex, increment)]

      ! get weights
      self%weights(1: self%subcells(i), i) = self%calculate_weights(&
              ! we cannot use the self%ids(1: self%subcells(i), i) directly as an indexer,
              ! as we need the one missing value (lower bound) for the first cell
              sourceCoord%values([(j, j = lowerSourceIndex + offset1, upperSourceIndex + offset2, increment)]), &
              lowerTargetBound, upperTargetBound&
      )
      if (any(le(self%weights(1: self%subcells(i), i), 0.0_dp))) then
        log_error(*) 'combine_weights: for targetCell:', i, &
                ' there are negative or zero-valued weights at ids: ', &
                self%ids(1: self%subcells(i), i)
        stop 1
      end if

      log_trace(*) 'compute_weights_1d_: for cell: ', i, ' with bounds: ', lowerTargetBound, upperTargetBound, &
        ' found ', self%subcells(i), ' subcells with ids: ', self%ids(1: self%subcells(i), i), &
        ' and weights: ', self%weights(1: self%subcells(i), i)

    end do

  end subroutine compute_weights_1d_

   function calculate_weights(self, sourceValues, lowerTargetBound, upperTargetBound) result(weights)
    !< this routine generates calculates the weights from the values

    !> the CoordUpscaler class
    class(CoordUpscaler), intent(in) :: self
    !> their ascending order is ensured
    real(dp), dimension(:), intent(in) :: sourceValues
    !> their ascending order is not ensured
    real(dp), intent(in) :: lowerTargetBound, upperTargetBound
    !> the weights slice for the sourceCell
    real(dp), dimension(size(sourceValues) - 1) :: weights

    integer(i4) :: nValues
    real(dp), dimension(size(sourceValues)) :: tempValues

    nValues = size(sourceValues)
    tempValues(:) = sourceValues(:)
    ! trim lower bound of first cell and upper bound of last cell to target cell bounds
    tempValues(1) = lowerTargetBound
    tempValues(nValues) = upperTargetBound
    ! get the weights by calculating the difference of the neighboring cells divided by total width
    weights(:) = abs(tempValues(2:nValues) - tempValues(1:nValues-1)) / (upperTargetBound - lowerTargetBound)

  end function calculate_weights
  
   function get_lower_source_index_asc_un(targetBound, vals, cachedIndex, step) result(index)
    real(dp), intent(in) :: targetBound
    real(dp), dimension(:), intent(in) :: vals
    integer(i8), intent(in) :: cachedIndex
    real(dp), intent(in), optional :: step
    integer(i8) :: index

    index = cachedIndex
    ! do while (vals(index + 1_i8) <= targetBound)
    do while ((targetBound - vals(index + 1_i8)) > (maxTolerance * (-1.0_dp)))
      index = index + 1_i8
    end do

  end function get_lower_source_index_asc_un

   function get_upper_source_index_asc_un(targetBound, vals, cachedIndex, step) result(index)
    real(dp), intent(in) :: targetBound
    real(dp), dimension(:), intent(in) :: vals
    integer(i8), intent(in) :: cachedIndex
    real(dp), intent(in), optional :: step
    integer(i8) :: index

    index = cachedIndex
    ! do while (vals(index + 1_i8) < targetBound)
    do while ((targetBound - vals(index + 1_i8)) > maxTolerance)
      index = index + 1_i8
    end do

  end function get_upper_source_index_asc_un

   function get_lower_source_index_desc_un(targetBound, vals, cachedIndex, step) result(index)
    real(dp), intent(in) :: targetBound
    real(dp), dimension(:), intent(in) :: vals
    integer(i8), intent(in) :: cachedIndex
    real(dp), intent(in), optional :: step
    integer(i8) :: index

    index = cachedIndex
    ! do while (vals(index) <= targetBound)
    do while ((targetBound - vals(index)) > (maxTolerance * (-1.0_dp)))
      index = index - 1_i8
    end do

  end function get_lower_source_index_desc_un

   function get_upper_source_index_desc_un(targetBound, vals, cachedIndex, step) result(index)
    real(dp), intent(in) :: targetBound
    real(dp), dimension(:), intent(in) :: vals
    integer(i8), intent(in) :: cachedIndex
    real(dp), intent(in), optional :: step
    integer(i8) :: index

    index = cachedIndex
    ! do while (vals(index) < targetBound)
    do while ((targetBound - vals(index)) > maxTolerance)
      index = index - 1_i8
    end do

  end function get_upper_source_index_desc_un

   function get_lower_source_index_asc_eq(targetBound, vals, cachedIndex, step) result(index)
    real(dp), intent(in) :: targetBound
    real(dp), dimension(:), intent(in) :: vals
    integer(i8), intent(in) :: cachedIndex
    real(dp), intent(in), optional :: step
    integer(i8) :: index, floored
    real(dp) :: advance

    ! floor means largest integer smaller than or equal to X - exactly, what we need
    advance = (targetBound - vals(cachedIndex)) / step
    floored = nint(advance, kind=i8)
    if (abs(real(floored, kind=dp) - advance) > maxTolerance ) then
      floored = floor(advance, kind=i8)
    end if
    index = cachedIndex + floored

  end function get_lower_source_index_asc_eq

   function get_upper_source_index_asc_eq(targetBound, vals, cachedIndex, step) result(index)
    real(dp), intent(in) :: targetBound
    real(dp), dimension(:), intent(in) :: vals
    integer(i8), intent(in) :: cachedIndex
    real(dp), intent(in), optional :: step
    integer(i8) :: index, floored
    real(dp) :: advance

    ! we need largest integer smaller than advance
    advance = (targetBound - vals(cachedIndex)) / step
    floored = nint(advance, kind=i8)
    if (abs(real(floored, kind=dp) - advance) > maxTolerance ) then
      index = cachedIndex + floor(advance, kind=i8)
    else
      index = cachedIndex + floored - 1_i8
    end if

  end function get_upper_source_index_asc_eq

   function get_lower_source_index_desc_eq(targetBound, vals, cachedIndex, step) result(index)
    real(dp), intent(in) :: targetBound
    real(dp), dimension(:), intent(in) :: vals
    integer(i8), intent(in) :: cachedIndex
    real(dp), intent(in), optional :: step
    integer(i8) :: index, floored
    real(dp) :: advance

    ! floor means largest integer smaller than or equal to X - exactly, what we need
    advance = (targetBound - vals(cachedIndex+1_i8)) / abs(step)
    floored = nint(advance, kind=i8)
    if (abs(real(floored, kind=dp) - advance) > maxTolerance ) then
      floored = floor(advance, kind=i8)
    end if
    index = cachedIndex - floored

  end function get_lower_source_index_desc_eq

   function get_upper_source_index_desc_eq(targetBound, vals, cachedIndex, step) result(index)
    real(dp), intent(in) :: targetBound
    real(dp), dimension(:), intent(in) :: vals
    integer(i8), intent(in) :: cachedIndex
    real(dp), intent(in), optional :: step
    integer(i8) :: index, floored
    real(dp) :: advance

    ! we need largest integer smaller than advance
    advance = (targetBound - vals(cachedIndex+1_i8)) / abs(step)
    floored = nint(advance, kind=i8)
    if (abs(real(floored, kind=dp) - advance) > maxTolerance ) then
      index = cachedIndex - floor(advance, kind=i8)
    else
      index = cachedIndex - floored + 1_i8
    end if

  end function get_upper_source_index_desc_eq

  subroutine combine_weights(self, coordIds, coordNSourceCells)
    class(CoordUpscaler), intent(inout) :: self
    integer(i4), dimension(:), intent(in) :: coordIds
    integer(i8), dimension(:), intent(in) :: coordNSourceCells

    type(CoordUpscalerPointer), dimension(size(coordIds)) :: upscalers
    integer(i8) :: nTargetCells, iCell, iSubcell, idMultiplier, nSourceCells, prodMaxSubcells
    integer(i4) :: iId, nIds
    integer(i8), dimension(size(coordIds, kind=i8)) :: coordTargetCells, coordTargetIdHelper, coordTargetIds, dimNSubcells, &
            dimSubcellIds, dimSubcellIdHelper, dimSourceIdHelper, maxSubcells
    integer(i8), dimension(:), allocatable :: tempIds, idAddon, subcellIds
    real(dp), dimension(:), allocatable :: tempWeights

    log_debug(*) 'combine_weights: combining weights of Upscalers:'
    ! loop over all coordIds and set the pointers to the CoordUpscalers
    nIds = size(coordIds)
    do iId=1, nIds
      log_debug(*) 'combine_weights:     ', iId, ': ', trim(MPR_COORD_UPSCALERS(coordIds(iId))%name)
      upscalers(iId)%coord_p => MPR_COORD_UPSCALERS(coordIds(iId))
    end do

    ! infer our number of target cells for each coordinate
    coordTargetCells = [(upscalers(iId)%coord_p%get_n_target_cells(), iId=1, nIds)]
    ! infer the number of target cells for the target grid
    nTargetCells = product(coordTargetCells)
    ! infer our number of maximum subcells needed
    maxSubcells = [(upscalers(iId)%coord_p%get_max_subcells(), iId=1, nIds)]
    prodMaxSubcells = product(maxSubcells)

    ! initialize our target arrays
    log_trace(*) 'combine_weights: init_values with maxSubcells ', prodMaxSubcells, ' and nTargetCells ', nTargetCells
    call self%init_values(prodMaxSubcells, nTargetCells)

    if (all([(.not. upscalers(iId)%coord_p%doNeedWeights, iId=1, nIds)])) then
      ! fast calculation for special case with equidistant cells
      self%doNeedWeights = .false.
      ! self%weights(:, :) = product([(upscalers(iId)%coord_p%weights(1_i8, 1_i8), iId=1, nIds)])
      ! ! loop over all the cells of the target
      ! nTotalSubcells(:) = maxSubcells * coordTargetCells
      ! sliceMins = []
      ! allocate(sliceMins(, product([(upscalers(iId)%coord_p%subcells(1_i8), iId=2, nIds)])))
      ! addOns = []
      ! sliceMaxs =[]
      !   do iBlock=1, product(maxSubcells(2:nIds))
      ! do iId=1, nIds
      !     iMax = iMin + maxSubcells(1_i4) - 1_i4
      !     iMin = (iBlock - 1_i4) * maxSubcells(1_i4) + 1_i4
      !             upscalers(iId)%coord_p%ids(:, :) + addOn(iSlice, iId)
      !     self%ids(iMin: iMax, sliceMins(iSlice, iId): sliceMaxs(iSlice, iId)) = &
      !   end do
      ! end do
      ! self%subcells(:) = prodMaxSubcells
    end if
    ! initialize the helper arrays
    allocate(tempWeights(prodMaxSubcells))
    allocate(tempIds(prodMaxSubcells))

    ! init the id helper (first value is always 1)
    dimSourceIdHelper(1) = 1_i8
    coordTargetIdHelper(1) = 1_i8
    do iId=2, nIds
      ! e.g. coordTargetCells = (4,16,16) -> coordTargetIds = (1,4,16)
      dimSourceIdHelper(iId) = product(coordNSourceCells(1:iId-1))
      coordTargetIdHelper(iId) = product(coordTargetCells(1:iId-1))
    end do

    ! loop over all the cells of the target
    do iCell=1, nTargetCells
      ! translate the iCell into a vector of Ids (size nId), getting the Id of each singular coordinate
      coordTargetIds = get_target_ids(coordTargetIdHelper, iCell)
      log_trace(*) 'combine_weights: targetCell:', iCell, ' with coordTargetIds: ', coordTargetIds
      ! get the number of subcells for each singular coordinate
      dimNSubcells = [(upscalers(iId)%coord_p%subcells(coordTargetIds(iId)), iId=1, nIds)]
      log_trace(*) 'combine_weights: targetCell:', iCell, ' with dimNSubcells: ', dimNSubcells
      ! get the total of resulting subcells for the target coordinate
      self%subcells(iCell) = product(dimNSubcells)
      if (self%subcells(iCell) > 0_i4) then
        ! get the index helper for the subcells
        dimSubcellIdHelper = 1_i8
        do iId=2, nIds
          ! e.g. coordTargetCells = (4,16,16) -> coordTargetIds = (1,4,16)
          dimSubcellIdHelper(iId) = product(dimNSubcells(1:iId-1))
        end do
        ! init the vectors of ids and weights for the current iCell
        tempWeights = 0.0_dp
        tempIds = 0_i8
        ! now loop over each subcell and get the correct target id and weight
        do iSubcell=1, self%subcells(iCell)
          ! translate the iSubcell into a vector of SubcellIds (size nId),
          ! getting the SubcellId of each singular coordinate
          dimSubcellIds = get_target_ids(dimSubcellIdHelper, iSubcell)
          tempWeights(iSubcell) = 1.0_dp
          do iId=1, nIds
            tempWeights(iSubcell) = tempWeights(iSubcell) * &
                    upscalers(iId)%coord_p%weights(dimSubcellIds(iId), coordTargetIds(iId))
            tempIds(iSubcell) = tempIds(iSubcell) + &
                    (upscalers(iId)%coord_p%ids(dimSubcellIds(iId), coordTargetIds(iId)) - 1_i8) * dimSourceIdHelper(iId)
            log_subtrace(*) 'combine_weights: for iID', iId, ' and upscaler ', &
                    upscalers(iId)%coord_p%ids(dimSubcellIds(iId), coordTargetIds(iId))
            log_subtrace(*) 'combine_weights: got id ', tempIds(iSubcell), ' and weight ', tempWeights(iSubcell)
          end do
          tempIds(iSubcell) = tempIds(iSubcell) + 1_i8
        end do
        log_trace(*) 'combine_weights: targetCell:', iCell, ' subcells', self%subcells(iCell)
        if (self%subcells(iCell) > 3_i8) then
          log_trace(*) 'combine_weights: targetCell:', iCell, ' weights', tempWeights(1:3), '...'
          log_trace(*) 'combine_weights: targetCell:', iCell, ' ids', tempIds(1:3), '...'
        else
          log_trace(*) 'combine_weights: targetCell:', iCell, ' weights', tempWeights(1:self%subcells(iCell))
          log_trace(*) 'combine_weights: targetCell:', iCell, ' ids', tempIds(1:self%subcells(iCell))
        end if
        if (any(le(tempWeights(1:self%subcells(iCell)), 0.0_dp))) then
          log_error(*) 'combine_weights: for targetCell:', iCell, &
                  ' there are negative or zero-valued weights at ids: ', &
                  tempIds(1:self%subcells(iCell))
          stop 1
        end if
        self%weights(1:self%subcells(iCell), iCell) = tempWeights(1:self%subcells(iCell))
        self%ids(1:self%subcells(iCell), iCell) = tempIds(1:self%subcells(iCell))
      end if
    end do

    deallocate(tempWeights, tempIds)

    log_debug(*) 'combine_weights: flag doNeedWeights: ', self%doNeedWeights

    ! set the mapMethod
    self%mapMethod = upscalers(1)%coord_p%mapMethod
    do iId=2, nIds
      if (self%mapMethod /= upscalers(iId)%coord_p%mapMethod) then
        self%mapMethod = 'multiple mapping methods for each coordinate'
        exit
      end if
    end do

  end subroutine combine_weights

  subroutine from_other(self, other)
    class(CoordUpscaler), intent(inout) :: self
    class(CoordUpscaler), intent(in) :: other

    self%subcells = other%get_subcells()
    self%ids = other%get_ids()
    self%weights = other%get_weights()
    self%mapMethod = other%mapMethod

  end subroutine from_other

  function get_target_ids(coordTargetIdHelper, iCell) result(coordTargetIds)
    !< compute the ids in each subdimension for a given id in the target coordinate
    integer(i8), dimension(:), intent(in):: coordTargetIdHelper
    integer(i8), intent(in ):: iCell
    integer(i8), dimension(size(coordTargetIdHelper)):: coordTargetIds

    integer(i8) :: iId, remainder

    remainder = iCell
    do iId=size(coordTargetIdHelper), 1, -1
      coordTargetIds(iId) = int((remainder-1_i8)/coordTargetIdHelper(iId)) + 1_i8
      remainder = mod((remainder-1_i8), coordTargetIdHelper(iId)) + 1_i8
    end do


  end function get_target_ids

  subroutine read_weights(self, fromWeightFile, subcellIdsFieldName, weightsFieldName, nSubcellsFieldName)
    !< read weights from netcdf file
    !< either use a SCRIP-based format or a custom format with additional nsubsell information
    class(CoordUpscaler), intent(inout) :: self
    character(*), intent(in) :: fromWeightFile
    character(*), intent(in), optional :: subcellIdsFieldName
    character(*), intent(in), optional :: weightsFieldName
    character(*), intent(in), optional :: nSubcellsFieldName

    type(NcDataset) :: nc
    type(NcVariable) :: ncVar
    integer(i4), dimension(:), allocatable :: varShape
    integer(i8), dimension(:), allocatable :: grid2_add_map1, grid1_add_map1
    real(dp), dimension(:), allocatable :: temp_dp
    real(dp), dimension(:,:), allocatable :: wts_map1
    integer(i4), dimension(:), allocatable :: temp_i4
    integer(i8) :: tempIndexer, iSubcell
    integer(i8) :: initSize

    log_debug(*) "read_weights: Initializing CoordUpscaler '", trim(self%name), "' from file."

    ! read Dataset and store the variable of interest
    nc = NcDataset(fromWeightFile, "r")
    ! infer the data format
    if (nc%hasVariable(trim(scripWeightsName)) .and. &
            nc%hasVariable(trim(scripSrcAddressName)) .and. &
            nc%hasVariable(trim(scripDstAddressName))) then
      ! SCRIP format
      ncVar = nc%getVariable(trim(scripWeightsName))
      call ncVar%getData(wts_map1)
      ncVar = nc%getVariable(trim(scripSrcAddressName))
      call ncVar%getData(grid1_add_map1)
      ncVar = nc%getVariable(trim(scripDstAddressName))
      call ncVar%getData(grid2_add_map1)
      ! read data from netcdf file, now init the value to the square root of the length of the sparse vector
      initSize = int(sqrt(real(size(grid2_add_map1, kind=i8))), kind=i8)
      call self%init_values(initSize, initSize)

      call convert_weights_format(wts_map1, grid1_add_map1, grid2_add_map1, self%weights, self%ids, self%subcells)
    else if (.not. present(nSubcellsFieldName) .and. &
      .not. present(subcellIdsFieldName) .and. &
      .not. present(weightsFieldName)) then
      log_error(*) 'read_weights: Neither recognized valid SCRIP format ', &
      '(with variables ', trim(scripSrcAddressName), ', ', trim(scripDstAddressName), &
              ', ', trim(scripWeightsName), ')', &
      'nor supply valid arguments for nSubcellsFieldName, subcellIdsFieldName and weightsFieldName.'
      call nc%close()
      stop 1
    else
      if (nc%hasVariable(trim(nSubcellsFieldName))) then
        ! the n_subcells vector is found in the nc file
        ncVar = nc%getVariable(trim(nSubcellsFieldName))
        varShape = ncVar%getShape()
        allocate(self%subcells(varShape(1)))
        call ncVar%getData(self%subcells)

        ! init the other target arrays
        allocate(self%ids(maxval(self%subcells), size(self%subcells, kind=i8)))
        allocate(self%weights(maxval(self%subcells), size(self%subcells, kind=i8)))
        self%ids = 0_i4
        self%weights = 0.0_dp

      else
        log_error(*) 'read_weights: ', trim(nSubcellsFieldName), &
                " is not a variable in the nc file ", trim(fromWeightFile)
        call nc%close()
        stop 1
      end if
      if (nc%hasVariable(trim(subcellIdsFieldName))) then
        ! the coordinate is actually a variable in the NcDataset
        ncVar = nc%getVariable(trim(subcellIdsFieldName))
        varShape = ncVar%getShape()
        allocate(temp_i4(varShape(1)))
        call ncVar%getData(temp_i4)
        tempIndexer = 1_i8
        do iSubcell=1_i8, size(self%subcells, kind=i8)
            self%ids(1:self%subcells(iSubcell), iSubcell) = temp_i4(&
                    tempIndexer: tempIndexer + self%subcells(iSubcell)-1)
          tempIndexer = tempIndexer + self%subcells(iSubcell)
        end do
        deallocate(temp_i4)
      else
        log_error(*) trim(subcellIdsFieldName), " is not a variable in the nc file ", trim(fromWeightFile)
          call nc%close()
        stop 1
      end if
      if (nc%hasVariable(trim(weightsFieldName))) then
        ! the coordinate is actually a variable in the NcDataset
        ncVar = nc%getVariable(trim(weightsFieldName))
        varShape = ncVar%getShape()
        allocate(temp_dp(varShape(1)))
        call ncVar%getData(temp_dp)
        tempIndexer = 1_i8
        do iSubcell=1_i8, size(self%subcells, kind=i8)
            self%weights(1:self%subcells(iSubcell), iSubcell) = temp_dp(&
                    tempIndexer: tempIndexer + self%subcells(iSubcell)-1)
          if (any(le(self%weights(1:self%subcells(iSubcell), iSubcell), 0.0_dp))) then
            log_error(*) 'read_weights: for targetCell', iSubcell, &
                    ' there are negative or zero-valued weights at ids: ', &
                    self%ids(1:self%subcells(iSubcell), iSubcell)
            stop 1
          end if
          tempIndexer = tempIndexer + self%subcells(iSubcell)
        end do
        deallocate(temp_dp)
      else
        log_error(*) trim(weightsFieldName), " is not a variable in the nc file ", trim(fromWeightFile)
          call nc%close()
        stop 1
      end if
    end if

    ! read the map_method from the attribute
    if (nc%hasAttribute('map_method')) then
      call nc%getAttribute('map_method', self%mapMethod)
    else
      log_warn(*) 'read_weights: Did not find attribute "map_method" in file ', trim(fromWeightFile)
    end if
    call nc%close()

  end subroutine read_weights

  subroutine write_weights(self)
    !< write weights to netcdf file
    class(CoordUpscaler), intent(inout) :: self

    type(NcDataset) :: nc
    type(NcVariable) :: ncVar
    type(NcDimension) :: NcDimWeights, NcDimLinks
    character(maxStringLength) :: toWeightFile
    logical :: doesFileExist
    integer(i4), parameter :: numWeights = 1_i4
    integer(i8) :: numLinks, iCell, iLinks, indexMin, indexMax, i
    integer(i4) :: iChar
    integer(i8), dimension(:), allocatable :: srcAddress, dstAddress
    real(dp), dimension(:, :), allocatable :: remapMatrix
    character(8) :: date
    character(3) :: integer_string
    character(maxNameLength) :: dimName
    !character(10) :: time
    !character(5)  :: zone

    log_debug(*) "write_weights: Writing CoordUpscaler '", trim(self%name), "' to file."

    ! read Dataset and store the variable of interest
    toWeightFile = 'weights_'//trim(self%alias)//'.nc'
    ! check if file already exists, then pass
    !inquire(file=toWeightFile, exist=doesFileExist)
    doesFileExist = .false.
    if (.not. doesFileExist) then
      ! init the arrays
      numLinks = sum(self%subcells)
      ! open file
      nc = NcDataset(toWeightFile, "w")

      ! write dimensions to nc file
      if (numLinks <= int(huge(0_i4), i8)) then
        allocate(srcAddress(numLinks), dstAddress(numLinks), remapMatrix(numWeights, numLinks))
        iLinks = 1_i8
        ! build SCRIP sparse matrix
        do iCell=1, size(self%subcells, kind=i8)
          srcAddress(iLinks:iLinks+self%subcells(iCell)-1) = self%ids(1:self%subcells(iCell), iCell)
          dstAddress(iLinks:iLinks+self%subcells(iCell)-1) = iCell
          remapMatrix(1, iLinks:iLinks+self%subcells(iCell)-1) = self%weights(1:self%subcells(iCell), iCell)
          iLinks = iLinks+self%subcells(iCell)
        end do
        log_debug(*) "Write dimension '", trim(scripWeightsCoord), "' to file ", trim(nc%fname)
        NcDimWeights = nc%setDimension(trim(scripWeightsCoord), numWeights)
        ! write the integer dimensions
        log_debug(*) "Write dimension '", trim(scripLinksCoord), "' to file ", trim(nc%fname)
        NcDimLinks = nc%setDimension(trim(scripLinksCoord), int(numLinks, i4))
        ! now set everything related to the data
        log_debug(*) "Write variable '", trim(scripSrcAddressName), "' to file ", trim(nc%fname)
        ncVar = nc%setVariable(trim(scripSrcAddressName), 'i32', [NcDimLinks])
        call ncVar%setData(int(srcAddress, i4))
        log_debug(*) "Write variable '", trim(scripDstAddressName), "' to file ", trim(nc%fname)
        ncVar = nc%setVariable(trim(scripDstAddressName), 'i32', [NcDimLinks])
        call ncVar%setData(int(dstAddress, i4))
        log_debug(*) "Write variable '", trim(scripWeightsName), "' to file ", trim(nc%fname)
        ncVar = nc%setVariable(trim(scripWeightsName), "f64", [NcDimWeights, NcDimLinks])
        call ncVar%setData(remapMatrix)
        deallocate(srcAddress, dstAddress, remapMatrix)
        ! set the attributes
        call nc%setAttribute("title", trim(self%alias))
        call nc%setAttribute("normalization", 'fracarea')
        call nc%setAttribute("map_method", self%mapMethod)
        call date_and_time(date=date)
        call nc%setAttribute("history", 'Created: '//date)
        iChar = index(self%alias, '__to__')
        call nc%setAttribute("source_grid", self%alias(1:iChar-1_i4))
        call nc%setAttribute("dest_grid", trim(self%alias(iChar+6_i4:)))
        call nc%setAttribute("conventions", 'SCRIP')
        call nc%close()
      else
        ! set unlimited dimension
        log_warn(*) 'write_weights: cannot write Coordinate with length >', int(huge(0_i4), i8), ' to netcdf file.'
        ! NcDimLinks = nc%setDimension(trim(scripLinksCoord))
      end if

    end if


  end subroutine write_weights

  subroutine init_values(self, maxSubcells, nTargetCells)
    !< allocate and initialize the core arrays containing the information on weights and subcellIds
    class(CoordUpscaler), intent(inout) :: self
    integer(i8), intent(in) :: maxSubcells
    integer(i8), intent(in) :: nTargetCells

    allocate(self%subcells(nTargetCells))
    allocate(self%ids(maxSubcells, nTargetCells))
    allocate(self%weights(maxSubcells, nTargetCells))
    self%subcells = 0_i8
    self%ids = 0_i8
    self%weights = 0.0_dp

  end subroutine init_values

  function get_n_target_cells(self) result(nTargetCells)
    class(CoordUpscaler), intent(in) :: self
    integer(i8) :: nTargetCells

    if (allocated(self%subcells)) then
      nTargetCells = size(self%subcells, kind=i8)
    else
      nTargetCells = -1_i8
    end if
  end function get_n_target_cells

  function get_max_subcells(self) result(maxSubcells)
    class(CoordUpscaler), intent(in) :: self
    integer(i8) :: maxSubcells

    if (allocated(self%subcells)) then
      maxSubcells = maxval(self%subcells)
    else
      maxSubcells = -1_i8
    end if
  end function get_max_subcells

  function get_subcells(self) result(subcells)
    class(CoordUpscaler), intent(in) :: self
    integer(i8), dimension(:), allocatable :: subcells

    if (allocated(self%subcells)) then
      subcells = self%subcells
    end if
  end function get_subcells

  function get_ids(self) result(ids)
    class(CoordUpscaler), intent(in) :: self
    integer(i8), dimension(:,:), allocatable :: ids

    if (allocated(self%ids)) then
      ids = self%ids
    end if
  end function get_ids

  function get_weights(self) result(weights)
    class(CoordUpscaler), intent(in) :: self
    real(dp), dimension(:,:), allocatable :: weights

    if (allocated(self%weights)) then
      weights = self%weights
    end if
  end function get_weights

  function get_index_in_coord_upscaler(itemName, globalvaluesArg, aliasName) result(iCoord)
    !< returns the index of the CoordUpscaler of globalvaluesArg if its name or alias matches itemName
    !< optionally, a match against aliasName is also checked

    !> the CoordUpscaler name to be found
    character(*), intent(in), optional :: itemName
    !> the array of CoordUpscaler to be used, falls back to MPR_COORD_UPSCALERS if not provided
    type(CoordUpscaler), dimension(:), intent(in), optional, target :: globalvaluesArg
    !> alias name to be checked with the reference in addition to itemName
    character(*), intent(in), optional :: aliasName
    integer(i4) :: iCoord

    type(CoordUpscaler), dimension(:), pointer :: globalVector => null()

    integer(i4) :: nCoords

    if (present(globalvaluesArg)) then
      globalVector => globalvaluesArg
    else
      globalVector => MPR_COORD_UPSCALERS
    end if

    ! get number of initialized coordinates
    nCoords = get_n_initialized(globalVector)
    ! get the index
    do iCoord = 1, nCoords
      ! check for name and alias
      if (trim(itemName) == trim(globalVector(iCoord)%name)) then
        return
      else if (trim(itemName) == trim(globalVector(iCoord)%alias)) then
        return
      end if
      if (present(aliasName)) then
        if (trim(aliasName) == trim(globalVector(iCoord)%name)) then
          return
        else if (trim(aliasName) == trim(globalVector(iCoord)%alias)) then
          return
        end if
      end if
    end do

  end function get_index_in_coord_upscaler

  subroutine get_stats(self)
    !< prints some information on the CoordUpscaler
    !< useful for debugging or printing before raising error messages

    !> CoordUpscaler type-bound procedure
    class(CoordUpscaler), intent(in) :: self

    log_info("(1X,A,A,A,A,I0,A,A,A)") 'CoordUpscaler "', trim(self%name), '" with id ', self%id, &
            ' and alias "', trim(self%alias), '"'
    log_info(*) 'The flag "doNeedWeights" is: ', self%doNeedWeights

    if (allocated(self%subcells)) then
      log_info("(1X,A,A,I0,A,I0,A,I0)") 'There are ', size(self%subcells, kind=i8), &
              ' target cells with the number of subcells ranging from ', &
              minval(self%subcells), ' to ', maxval(self%subcells)
      log_info("(1X,A,I0,A)") size(pack(self%subcells, self%subcells==0), kind=i8), ' cells have zero subcells'
      log_info("(1X,A,A,I0,A,I0,A,I0)") 'There are ', size(self%weights,2, kind=i8), &
              ' target cells with a maximum of ', size(self%weights,1, kind=i8), ' subcells each in the weights array'
      log_info("(1X,A,A,F0.8,A,F0.8)") 'The valid weights are ranging from ', &
              minval(self%weights, self%weights>0.0_dp), ' to ', maxval(self%weights, self%weights>0.0_dp)
      log_info("(1X,A,A,I0,A,I0,A,I0)") 'There are ', size(self%ids,2, kind=i8), ' target cells with a maximum of ', &
              size(self%ids,1, kind=i8), ' subcells each in the weights array'
      log_info("(1X,A,A,I0,A,I0)") 'The valid ids are ranging from ', &
              minval(self%ids, self%ids>0_i4), ' to ', maxval(self%ids, self%ids>0_i4)
    else
      log_info(*) 'The values are not yet initialized.'
    end if


  end subroutine get_stats


  subroutine convert_weights_format(weightsIn, sourceId, targetId, weights, ids, subcells)
    !< translate the weights and links that are stored by SCRIP into
    !< data format needed by MPR

    real(dp), intent(in), dimension(:,:) :: weightsIn
    integer(i8), intent(in), dimension(:) :: sourceId
    integer(i8), intent(in), dimension(:) :: targetId
    real(dp), intent(inout), dimension(:,:), allocatable :: weights
    integer(i8), intent(inout), dimension(:,:), allocatable :: ids
    integer(i8), intent(inout), dimension(:), allocatable :: subcells

    integer(i8) :: iLink, lastLink, iSubcell
    real(dp) :: extremeWeight

    log_debug(*) 'convert_weights_format: Translating the weights and links/ids from SCRIP format to MPR format'
    ! init the counter of the current subcell index
    iSubcell = 0_i8
    ! init the lastLink to the first address of the target grid
    lastLink = targetId(1)
    ! loop over all the links that have weights
    ! TODO: think about fast track for grids with all cells having the same number of corners
    log_trace(*) 'sizes of input: sourceId (', shape(sourceId, kind=i8), '), targetId (', size(targetId, kind=i8), &
    ') and weightsIn', shape(weightsIn, kind=i8), ')'
    log_trace(*) 'sizes of output: ids (', shape(ids, kind=i8), '), subcells (', size(subcells, kind=i8), &
    ') and weights', shape(weights, kind=i8), ')'
    log_subtrace(*) 'sourceId:', sourceId
    log_subtrace(*) 'targetId:', targetId
    log_subtrace(*) 'weightsIn:', weightsIn

    do iLink=1, size(targetId)
      if (targetId(iLink) > size(subcells, kind=i8)) then
        log_debug(*) 'convert_weights_format: increasing arrays to accomodate more target cells from ', &
            size(weights, 2, kind=i8), ' to ', size(weights, 2, kind=i8)*2
        call resize_weights_array(weights, ids, subcells, &
                size(weights, 1, kind=i8), size(weights, 2, kind=i8)*2_i8)
      end if
      ! is there a new target grid address?
      if (targetId(iLink) /= lastLink) then
        ! update the number of subcells of the previous target address
        subcells(targetId(iLink-1)) = iSubcell
        ! continue counting the subcells from the last count
        iSubcell = subcells(targetId(iLink))
        ! update the lastlink
        lastLink = targetId(iLink)
      end if
      iSubcell = iSubcell + 1_i8

      if (iSubcell > size(weights, 1, kind=i8)) then
        log_debug(*) 'convert_weights_format: increasing arrays to accomodate more subcells from ', &
            size(weights, 1, kind=i8), ' to ', size(weights, 1, kind=i8)*2
        call resize_weights_array(weights, ids, subcells, &
                size(weights,1)*2_i8, size(weights,2, kind=i8))
      end if
      if (targetId(iLink) > size(weights, 2, kind=i8)) then
        log_debug(*) 'convert_weights_format: increasing arrays to accomodate more target cells from ', &
            size(weights, 2, kind=i8), ' to ', size(weights, 2, kind=i8)*2
        call resize_weights_array(weights, ids, subcells, &
                size(weights, 1, kind=i8), size(weights, 2, kind=i8)*2_i8)
      end if


      ! use the first order weights only
      weights(iSubcell, targetId(iLink)) = weightsIn(1, iLink)
      ! use the source address
      ids(iSubcell, targetId(iLink)) = sourceId(iLink)
    end do

    if (iSubcell /= 0_i8) then
      ! update the number of subcells of the previous target address
      subcells(targetId(iLink-1)) = iSubcell
    end if

    log_debug(*) 'convert_weights_format: Done translating weights and links'
    log_trace(*) 'convert_weights_format: Global sum of weights: ', sum(weightsIn(1, :))
    extremeWeight = maxval(weightsIn(1, :))
    log_trace(*) 'convert_weights_format: Global max of weights: ', extremeWeight
    if (extremeWeight > 1.0_dp) then
      log_error(*) 'SCRIP calculated weights > 1: ', extremeWeight, &
      '. Check if you provided polygons with corners not in a counter-clockwise order.'
      !' There is also a known bug in SCRIP with custom-made polygons.'
      stop 1
    end if
    extremeWeight = minval(weightsIn(1, :))
    log_trace(*) 'convert_weights_format: Global min of weights: ', extremeWeight
    if (extremeWeight < 0.0_dp) then
      log_error(*) 'SCRIP calculated weights < 0: ', extremeWeight, &
      '. Check if you provided polygons with corners not in a counter-clockwise order.'
      !' There is also a known bug in SCRIP with custom-made polygons.'
      stop 1
    end if

    call resize_weights_array(weights, ids, subcells, nCellsArg=count(subcells/=0_i8, kind=i8))

  end subroutine convert_weights_format

  subroutine resize_weights_array(weights, ids, subcells, maxSubcellsArg, nCellsArg)
    !< change the size of the arrays; if increased, then fill with default values
    !< from weights(m,n), ids(m,n), subcells(n) to
    !< case no args: m = maxval(subcells), n = size(subcells)
    !< case maxSubcellsArg: m = maxSubcellsArg, n = size(subcells)
    !< case nCellsArg: m = maxval(subcells), n = nCellsArg
    !< case maxSubcellsArg, nCellsArg: m = maxSubcellsArg, n = nCellsArg

    !> the array of weights (subcells, target cell)
    real(dp), intent(inout), dimension(:,:), allocatable :: weights
    !> the array of subcell ids (subcell ids, target cell)
    integer(i8), intent(inout), dimension(:,:), allocatable :: ids
    !> the array of the number of subcell ids (target cell)
    integer(i8), intent(inout), dimension(:), allocatable :: subcells
    !> the optional number of subcell ids to allocate space for
    integer(i8), intent(in), optional :: maxSubcellsArg
    !> the optional number of subcell ids to allocate space for
    integer(i8), intent(in), optional :: nCellsArg

    integer(i8) :: maxSubcells, nCells, iCoord1, iCoord2
    real(dp), dimension(:,:), allocatable :: tempWeights
    integer(i8), dimension(:,:), allocatable :: tempIds
    integer(i8), dimension(:), allocatable :: temp1D

    log_debug(*) 'resize_weights_array: Adapting the size of the weights, ids and subcells arrays'

    ! alternative, given by user, can also be higher than maxSubcells
    if (present(maxSubcellsArg)) then
      maxSubcells = maxSubcellsArg
    else
      ! the maximum number of subcells
      maxSubcells = maxval(subcells)
    end if

    if (present(nCellsArg)) then
      nCells = nCellsArg
    else
      nCells = size(subcells, kind=i8)
    end if

    if (maxSubcells > size(weights, 1, kind=i8)) then
      ! if the array is enlarged, index all existing values
      iCoord1 = size(weights, 1, kind=i8)
    else
      ! if the array is made smaller, index values up to new size
      iCoord1 = maxSubcells
    end if
    if (nCells > size(weights, 2, kind=i8)) then
      iCoord2 = size(weights, 2, kind=i8)
    else
      iCoord2 = nCells
    end if

    ! working on the weights array
    if (maxSubcells /= size(weights, 1, kind=i8) .or. nCells /= size(weights, 2, kind=i8)) then
      allocate(tempWeights(maxSubcells, nCells))
      tempWeights(1:iCoord1, 1:iCoord2) = weights(1:iCoord1, 1:iCoord2)
      log_trace(*) 'resizing the weights array from (', shape(weights, kind=i8), &
              ') to size (', shape(tempWeights, kind=i8), ')'
      if (maxSubcells > size(weights, 1, kind=i8)) then
        tempWeights(iCoord1+1:, :) = 0_i8
      end if
      if (nCells > size(weights, 2, kind=i8)) then
        tempWeights(:, iCoord2+1:) = 0_i8
      end if
      call move_alloc(tempWeights, weights)
    end if

    ! working on the ids array
    if (maxSubcells /= size(ids, 1, kind=i8) .or. nCells /= size(ids, 2, kind=i8)) then
      allocate(tempIds(maxSubcells, nCells))
      tempIds(1:iCoord1, 1:iCoord2) = ids(1:iCoord1, 1:iCoord2)
      log_trace(*) 'resizing the ids array from (', shape(ids, kind=i8), &
              ') to size (', shape(tempIds, kind=i8), ')'
      if (maxSubcells > size(ids, 1, kind=i8)) then
        tempIds(iCoord1+1:, :) = 0_i8
      end if
      if (nCells > size(ids, 2, kind=i8)) then
        tempIds(:, iCoord2+1:) = 0_i8
      end if
      call move_alloc(tempIds, ids)
    end if

    ! working on the subcells vector
    if (nCells /= size(subcells, kind=i8)) then
      allocate(temp1D(nCells))
      temp1D(1:iCoord2) = subcells(1:iCoord2)
      log_trace(*) 'resizing the subcells array from (', shape(subcells, kind=i8), &
              ') to size (', shape(temp1D, kind=i8)
      if (nCells > size(subcells, kind=i8)) then
        temp1D(iCoord2:) = 0_i8
      end if
      call move_alloc(temp1D, subcells)
    end if

  end subroutine resize_weights_array

end module mo_mpr_data_array_upscale

