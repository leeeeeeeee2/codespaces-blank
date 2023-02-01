!> \file mo_mpr_data_array.f90

!> \brief MPR related types

!> \details A wrapper around the NetCDF Fortran 90 interface.
!
!> \authors Robert Schweppe, Stephan Thober
!> \date Jan 2018
#include "flogging.h"
module mo_mpr_data_array

  use mo_kind, only : dp, i4, i8
  use mo_mpr_constants, only : defaultAlias, defaultCoordStagger, defaultDataArrayLimits, &
          defaultDataArrayToFile, maxStringLength, maxNoParameterPerDA, maxNoTemporalData, &
          upperBoundName, lowerBoundName, maxNameLength

  use mo_mpr_utils, only : get_index_in_vector, get_n_initialized, &
          get_mask_and_packed_data1D, get_mask_and_packed_data2D, &
          get_mask_and_packed_data3D, get_mask_and_packed_data4D
  use mo_mpr_constants, only : maxNoDataArrays, maxStringLength, maxNameLength
  use mo_mpr_data_array_upscale, only : UpscaleHelper, add_upscaler, MPR_UPSCALERS
  use mo_mpr_global_variables, only : CHECK_FOR_NODATAVALUE

  use mo_string_utils, only : num2str, Replace_Text, compress, index_word, replace_word
  use mo_netcdf, only : NcDataset, NcDimension, NcVariable
  use mo_constants, only : nodata_dp
  use mo_mpr_input_field_container, only : InputFieldContainer
  use mo_mpr_transfer_func, only : mprTransferFunctionTable
  use mo_utils, only : eq, ne
  use mo_mpr_coordinate, only : Coordinate, MPR_COORDINATES, &
          CoordinatePointer, add_coordinate, get_index_in_coordinate, create_upscaler_name
  use mo_mpr_util_types, only : MprBaseType
  use mo_mpr_file, only : filenameNamelistMprDefault
  use mo_mpr_parameters, only : MPR_PARAMETERS
  use mo_orderpack, only : sort_index
  use flogging
  use mo_mpr_transfer_func, only: transfer_function_1
  use mo_mpr_transfer_func, only: transfer_function_2

  implicit none

  private

  public :: add_data_array, delete_parents

  ! --------------------------------------------------------------------------------------

  type, public :: TransferHelper
    private
    character(maxNameLength) :: name
    integer(i4) :: id
    logical, public :: doInitFromFile
    type(CoordinatePointer), dimension(:), allocatable :: targetCoords
    character(maxStringLength) :: fileName
    character(maxNameLength), public, allocatable, dimension(:) :: inputFieldNames
    character(maxStringLength), public :: transferFuncName
    character(maxNameLength) :: transferFuncLabel
    real(dp), dimension(:), allocatable :: globalParameters

  contains
    private
    procedure, public :: init_and_execute_transfer
    procedure, public :: reset_transferHelper
    procedure :: check_and_set_input_field_names
    procedure :: set_input_field_containers
    procedure :: check_data_arrays
    procedure :: check_transfer_func_args
    procedure :: concatenate_input_fields
    procedure :: execute_transfer
    procedure :: read_nc
    procedure :: call_transfer_func
    procedure, public :: func_str_to_funcName
    procedure :: get_stats

  end type TransferHelper

  interface TransferHelper
    procedure newTransferHelper
  end interface TransferHelper

  type, public, extends(MprBaseType) :: DataArray
    private
    ! helper for executing the transfer function
    type(TransferHelper), public :: transferHelper

    ! helper for executing the upscaling
    character(maxNameLength), allocatable, dimension(:) :: upscaleOperatorNames
    type(UpscaleHelper), pointer :: upscaleHelper => null()

    ! other
    real(dp), dimension(2) :: limits
    character(maxNameLength), allocatable, dimension(:) :: deleteFieldNames

    ! finalize
    logical, public :: toFile

    ! main
    character(maxNameLength), allocatable, dimension(:) :: coordNames
    type(CoordinatePointer), allocatable, dimension(:), public :: coords
    real(dp), dimension(:), allocatable, public :: data
    logical, dimension(:), allocatable, public :: reshapedMask

  contains
    private
    procedure, public :: execute
    procedure, public :: write
    procedure, public :: is_finalized
    procedure, public :: reset
    procedure, public :: add_delete_items
    procedure :: get_data_shape

  end type DataArray

  interface DataArray
    procedure newDataArray
  end interface DataArray

  type(DataArray), dimension(maxNoDataArrays), public, target :: MPR_DATA_ARRAYS ! all data arrays

  type, public :: TemporalData
    private
    real(dp), dimension(:), allocatable, public :: data
  end type TemporalData

  type, public :: TemporalMask
    private
    logical, dimension(:), allocatable, public :: mask
  end type TemporalMask

  type(TemporalData), dimension(maxNoTemporalData), public, target :: MPR_TEMPORAL_DATA ! all temporal data
  type(TemporalMask), dimension(maxNoTemporalData), public, target :: MPR_TEMPORAL_MASK ! all temporal mask
  ! --------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------
contains

  function newDataArray(name, id, fromDataArraysArg, fromFile, transferFuncName, transferFuncLabel, fromParameters, &
          coordNames, upscaleOperatorNames, limitsArg, toFileArg, fromData, fromReshapedMask, fromCoords)
    character(*), intent(in) :: name
    integer(i4), intent(in) :: id
    character(*), intent(in), dimension(:), optional :: fromDataArraysArg
    character(*), intent(in), optional :: fromFile
    character(*), intent(in), optional :: transferFuncName
    character(*), intent(in), optional :: transferFuncLabel
    real(dp), intent(in), dimension(:), optional :: fromParameters
    character(*), intent(in), dimension(:), optional :: coordNames
    character(*), intent(in), dimension(:), optional :: upscaleOperatorNames
    real(dp), dimension(2), intent(in), optional :: limitsArg
    logical, intent(in), optional :: toFileArg
    real(dp), intent(in), dimension(:), optional :: fromData
    logical, intent(in), dimension(:), optional :: fromReshapedMask
    type(CoordinatePointer), intent(in), dimension(:), optional :: fromCoords
    type(DataArray) :: newDataArray


    integer(i4) :: iCoord, iError
    real(dp), dimension(2) :: limits = defaultDataArrayLimits
    logical :: toFile = defaultDataArrayToFile

    newDataArray%name = trim(name)
    newDataArray%id = id
    newDataArray%is_initialized = .true.
    allocate(newDataArray%deleteFieldNames(0))

    if (present(coordNames)) then
      if (size(coordNames) > 0) then
        newDataArray%coordNames = coordNames
        allocate(newDataArray%coords(size(coordNames)))
        do iCoord = 1, size(newDataArray%coordNames)
          call newDataArray%coords(iCoord)%set_coordinate_pointer(newDataArray%coordNames(iCoord), iError)
        end do
      end if
    end if
    if (present(upscaleOperatorNames)) then
      if (size(upscaleOperatorNames) > 0_i4) then
        newDataArray%upscaleOperatorNames = upscaleOperatorNames
      end if
    end if

    if (present(fromData) .and. present(fromReshapedMask) .and. present(fromCoords)) then
      ! directly initialized from data
      newDataArray%data = fromData
      newDataArray%reshapedMask = fromReshapedMask
      newDataArray%coords = fromCoords
    else
      ! the TransferHelper is initialized here
      newDataArray%transferHelper = TransferHelper(newDataArray%name, newDataArray%id, fromDataArraysArg, &
              fromFile, transferFuncName, transferFuncLabel, fromParameters)
    end if

    if (present(limitsArg)) then
      newDataArray%limits(:) = limitsArg(:)
    else
      newDataArray%limits(:) = limits(:)
    end if

    if (present(toFileArg)) then
      newDataArray%toFile = toFileArg
    else
      newDataArray%toFile = toFile
    end if
    log_subtrace(*) 'Finished initialization of newDataArray'

  end function newDataArray

  subroutine execute(self)
    class(DataArray), intent(inout) :: self

    type(CoordinatePointer), allocatable, dimension(:) :: sourceCoords
    type(CoordinatePointer), allocatable, dimension(:) :: coords
    real(dp), dimension(:), allocatable :: data
    logical, dimension(:), allocatable :: reshapedMask
    type(UpscaleHelper) :: upscaleHelper_
    integer(i4) :: id, iCoord, iError
    character(maxNameLength) :: upscalerName

    ! if already finalized skip
    if (self%is_finalized()) return

    ! initialize the data and get the results of the initializer
    call self%transferHelper%init_and_execute_transfer(data, reshapedMask, coords)

    ! only print status if used a proper function
    if (.not. any(reshapedMask)) then
      log_warn("(1X,A,A,A,A)") "Only missing values detected for DataArray '", trim(self%name), "'."
    else if (trim(self%transferHelper%transferFuncName) /= trim(defaultAlias)) then
      log_info("(1X,A,A,A,A,F0.6,A,F0.6)") "Just transferred DataArray '", trim(self%name), &
              "' with value range ", minval(data), " - ", maxval(data)
    end if
    ! correct limits
    if (ne(self%limits(1), nodata_dp)) then
      where(data < self%limits(1)) data = self%limits(1)
    end if
    if (ne(self%limits(2), nodata_dp)) then
      where(data > self%limits(2)) data = self%limits(2)
    end if

    ! set the things that are known only now that we know the source coords
    if (.not. allocated(self%coordNames)) then
      ! set the pointers to the dimensions
      self%coords = coords
      ! turn off the upscaling, create the upscaleHelper
      id = get_index_in_vector(trim(defaultAlias), MPR_UPSCALERS)
      upscaleHelper_ = UpscaleHelper(&
              name = trim(defaultAlias), &
              id = id, &
              doUpscale = .false.&
              )
    else
      ! create the source_dim pointers
      sourceCoords = coords
      ! check if all coords are associated
      do iCoord = 1, size(self%coordNames)
        if (.not. associated(self%coords(iCoord)%coord_p)) then
          ! if not associated (Coordinate was not yet known during init), try that now and raise if it fails
          call self%coords(iCoord)%set_coordinate_pointer(self%coordNames(iCoord), iError)
          if (iError == 1_i4) then
            log_error(*) 'execute: Could not set coordinate "', trim(self%coordNames(iCoord)) ,&
                    '" for DataArray "', trim(self%name), '"'
            stop 1
          end if
        end if
      end do

      ! create the upscaleHelper name
      upscalerName = create_upscaler_name(sourceCoords, self%coords, self%upscaleOperatorNames)

      ! create the upscaleHelper
      id = get_index_in_vector(upscalerName, MPR_UPSCALERS)
      upscaleHelper_ = UpscaleHelper(&
              name = upscalerName, &
              id = id, &
              doUpscale = .true., &
              targetCoords = self%coords, &
              sourceCoords = sourceCoords, &
              upscaleOperatorNames = self%upscaleOperatorNames &
              )
    end if

    ! add it to the upscaleHelper list
    call add_upscaler(upscaleHelper_, id)
    ! and link it
    self%upscaleHelper => MPR_UPSCALERS(id)

    ! do the upscaling
    if (self%upscaleHelper%doUpscale) call self%upscaleHelper%execute(data, reshapedMask)
    ! set as self
    call move_alloc(data, self%data)
    call move_alloc(reshapedMask, self%reshapedMask)

  end subroutine execute

  function is_finalized(self)
    class(DataArray), intent(in) :: self
    logical :: is_finalized

    is_finalized = allocated(self%data)
  end function is_finalized

  subroutine reset(self)
    class(DataArray), intent(inout) :: self

    if (self%is_initialized) then
      self%name = ''
      self%id=0
      self%is_initialized = .false.
      call self%transferHelper%reset_transferHelper()

      ! upscaling
      if (allocated(self%upscaleOperatorNames)) deallocate(self%upscaleOperatorNames)
      self%upscaleHelper => null()
      ! main
      if (allocated(self%coordNames)) deallocate(self%coordNames)
      if (allocated(self%coords)) deallocate(self%coords)
      if (allocated(self%data)) deallocate(self%data)
      if (allocated(self%reshapedMask)) deallocate(self%reshapedMask)
      if (allocated(self%deleteFieldNames)) deallocate(self%deleteFieldNames)
    end if

  end subroutine reset

  function get_data_shape(self) result(dataShape)
    class(DataArray), intent(in) :: self
    integer(i8), dimension(:), allocatable :: dataShape

    integer(i4) :: iCoord

    allocate(dataShape(size(self%coords)))
    do iCoord = 1, size(dataShape)
      dataShape(iCoord) = self%coords(iCoord)%coord_p%count
    end do

  end function get_data_shape

  subroutine write(self, nc)
    class(DataArray), intent(inout) :: self
    class(NcDataset), intent(inout) :: nc

    integer(i4) :: iCoord
    type(NcDimension), dimension(:), allocatable :: NcDims
    type(NcVariable) :: ncVar
    character(maxNameLength), dimension(:), allocatable :: coordNames
    integer(i8), dimension(:), allocatable :: upscaledShape
    logical, dimension(:), allocatable :: tempMask1d
    logical, dimension(:, :), allocatable :: tempMask2d
    logical, dimension(:, :, :), allocatable :: tempMask3d
    logical, dimension(:, :, :, :), allocatable :: tempMask4d
    logical, dimension(:, :, :, :, :), allocatable :: tempMask5d

    allocate(NcDims(size(self%coords)))
    allocate(coordNames(size(self%coords)))

    upscaledShape = self%get_data_shape()
    ! write the dimensions of the current data array, start with last
    do iCoord = size(self%coords), 1, -1
      ! get the name
      if (self%coords(iCoord)%coord_p%is_polygon()) then
        ! SCRIP format
        coordNames(iCoord) = 'grid_size'
      else
        coordNames(iCoord) = trim(self%coords(iCoord)%coord_p%name)
      end if
      ! is the dimension already in the nc file? (set by previous data array)
      if (nc%hasDimension(coordNames(iCoord))) then
        ! add it to our vector of ncCoordinates for the current array
        NcDims(iCoord) = nc%getDimension(coordNames(iCoord))
      else
        log_debug(*) "Write coordinate '", trim(coordNames(iCoord)), "' to file ", trim(nc%fname)
        call self%coords(iCoord)%coord_p%write_coordinate(nc, NcDims(iCoord))
        ! missing value (although there should be none)
        if (.not. self%coords(iCoord)%coord_p%is_polygon()) then
          ncVar = nc%getVariable(trim(coordNames(iCoord)))
          call ncVar%setAttribute("missing_value", nodata_dp)
        end if
      end if
    end do

    log_info(*) "Write variable '", trim(self%name), "' to file ", trim(nc%fname)
    log_debug(*) 'it has shape (', upscaledShape(:), ')'
    log_subtrace(*) 'write: it has the values ', self%data
    log_subtrace(*) 'write: it has the mask ', self%reshapedMask
    ! now set everything related to the data
    ncVar = nc%setVariable(trim(self%name), "f64", NcDims)
    call ncVar%setFillValue(nodata_dp)
    ! set the data
    ! the compiler wants it that way...
    ! this boilerplate with the temp_mask is required by the gfortran 8.1 compiler (and possibly onwards)
    select case(size(upscaledShape(:)))
    case(1)
      tempMask1d = reshape(self%reshapedMask, upscaledShape(1 : 1))
      call ncVar%setData(unpack(self%data, tempMask1d, nodata_dp))
      deallocate(tempMask1d)
    case(2)
      tempMask2d = reshape(self%reshapedMask, upscaledShape(1 : 2))
      call ncVar%setData(unpack(self%data, tempMask2d, nodata_dp))
      deallocate(tempMask2d)
    case(3)
      tempMask3d = reshape(self%reshapedMask, upscaledShape(1 : 3))
      call ncVar%setData(unpack(self%data, tempMask3d, nodata_dp))
      deallocate(tempMask3d)
    case(4)
      tempMask4d = reshape(self%reshapedMask, upscaledShape(1 : 4))
      call ncVar%setData(unpack(self%data, tempMask4d, nodata_dp))
      deallocate(tempMask4d)
    case(5)
      tempMask5d = reshape(self%reshapedMask, upscaledShape(1 : 5))
      call ncVar%setData(unpack(self%data, tempMask5d, nodata_dp))
      deallocate(tempMask5d)
    end select
    ! set the attributes
    call ncVar%setAttribute("long_name", trim(self%name))
    ! missing value
    call ncVar%setAttribute("missing_value", nodata_dp)


  end subroutine write

  subroutine add_delete_items(self, deleteFieldNames)
    class(DataArray), intent(inout) :: self
    character(maxNameLength), intent(in) :: deleteFieldNames(:)

    self%deleteFieldNames = deleteFieldNames

  end subroutine add_delete_items

  subroutine delete_parents(item)
    class(DataArray), intent(in) :: item
    character(maxNameLength), allocatable, dimension(:) :: deleteFieldNames
    integer(i4) :: j

    ! copy because array itself might be deleted
    deleteFieldNames = item%deleteFieldNames
    do j = 1, size(MPR_DATA_ARRAYS)
      if (any(deleteFieldNames == MPR_DATA_ARRAYS(j)%name)) then
        log_debug(*) 'delete_parents: deleting DataArray "', trim(MPR_DATA_ARRAYS(j)%name), '"'
        call MPR_DATA_ARRAYS(j)%reset()
      end if
    end do
    deallocate(deleteFieldNames)

  end subroutine delete_parents

  subroutine add_data_array(item, insertIndex)
    type(DataArray), intent(in) :: item
    integer(i4), intent(in) :: insertIndex

    ! check if there is still room in the vector for one more item
    if (insertIndex >= size(MPR_DATA_ARRAYS)) then
      log_error("(1X,A,A,I0,A,I0)") "Data array cannot be added at index ('", insertIndex, &
              "'). Vector of data arrays only has length: ", size(MPR_DATA_ARRAYS)
      stop 1
    end if

    ! add item
    MPR_DATA_ARRAYS(insertIndex) = item

  end subroutine add_data_array

  subroutine add_temporal_data(item)
    real(dp), dimension(:), allocatable, intent(inout) :: item
    integer(i4) :: i

    ! check if there is still room in the vector for one more item
    do i = 1, size(MPR_TEMPORAL_DATA)
      if (.not. allocated(MPR_TEMPORAL_DATA(i)%data)) then
        call move_alloc(item, MPR_TEMPORAL_DATA(i)%data)
        return
      end if
    end do
    log_error("(1X,A,A,I0,A)") "Temporal data cannot be added, there are already '", size(MPR_TEMPORAL_DATA), &
            "' temporal datasets contained"
    stop 1

  end subroutine add_temporal_data

  subroutine clear_temporal_data()
    integer(i4) :: i

    ! clear all the temporalData
    do i = 1, size(MPR_TEMPORAL_DATA)
      if (allocated(MPR_TEMPORAL_DATA(i)%data)) then
        deallocate(MPR_TEMPORAL_DATA(i)%data)
      else
        exit
      end if
    end do

  end subroutine clear_temporal_data

  function get_temporal_data_index() result(index)
    integer(i4) :: index
    integer(i4) :: i

    ! check if there is still room in the vector for one more item
    index = 0
    do i = 1, size(MPR_TEMPORAL_DATA)
      if (.not. allocated(MPR_TEMPORAL_DATA(i)%data)) then
        index = i - 1
        return
      end if
    end do
    index = 0

  end function get_temporal_data_index

  subroutine add_temporal_mask(item)
    logical, dimension(:), allocatable, intent(inout) :: item
    integer(i4) :: i

    ! check if there is still room in the vector for one more item
    do i = 1, size(MPR_TEMPORAL_MASK)
      if (.not. allocated(MPR_TEMPORAL_MASK(i)%mask)) then
        call move_alloc(item, MPR_TEMPORAL_MASK(i)%mask)
        return
      end if
    end do
    log_error("(1X,A,A,I0,A)") "Temporal mask cannot be added, there are already '", size(MPR_TEMPORAL_MASK), &
            "' temporal datasets contained"
    stop 1

  end subroutine add_temporal_mask

  subroutine clear_temporal_mask()
    integer(i4) :: i

    ! clear all the temporalData
    do i = 1, size(MPR_TEMPORAL_MASK)
      if (allocated(MPR_TEMPORAL_MASK(i)%mask)) then
        deallocate(MPR_TEMPORAL_MASK(i)%mask)
      else
        exit
      end if
    end do

  end subroutine clear_temporal_mask

  function get_temporal_mask_index() result(index)
    integer(i4) :: index
    integer(i4) :: i

    ! check if there is still room in the vector for one more item
    index = 0
    do i = 1, size(MPR_TEMPORAL_MASK)
      if (.not. allocated(MPR_TEMPORAL_MASK(i)%mask)) then
        index = i - 1
        return
      end if
    end do
    index = 0

  end function get_temporal_mask_index

  function newTransferHelper(name, id, fromDataArraysArg, fromFile, transferFuncName, transferFuncLabel, fromParameters)
    character(*), intent(in) :: name
    integer(i4), intent(in) :: id
    character(*), intent(in), dimension(:), optional :: fromDataArraysArg
    character(*), intent(in), optional :: fromFile
    character(*), intent(in), optional :: transferFuncName
    character(*), intent(in), optional :: transferFuncLabel
    real(dp), intent(in), dimension(:), optional :: fromParameters
    character(maxNameLength), dimension(:), allocatable :: fromDataArrays
    type(TransferHelper) :: newTransferHelper

    log_subtrace(*) 'Initialize type newTransferHelper'
    newTransferHelper%name = trim(name)
    newTransferHelper%id = id

    ! check for file read
    newTransferHelper%doInitFromFile = .false.
    if (present(fromFile)) then
      if (trim(fromFile) /= trim(defaultAlias)) then
        newTransferHelper%doInitFromFile = .true.
        newTransferHelper%fileName = trim(fromFile)

        if (present(fromDataArraysArg)) then
          if (size(fromDataArraysArg) > 0) then
            log_debug(*) 'Data array "', trim(newTransferHelper%name), '" uses "fromFile" and', &
                  ' "fromDataArrays": ', fromDataArraysArg(:)
            call newTransferHelper%check_and_set_input_field_names(fromDataArraysArg)
            allocate(fromDataArrays(size(fromDataArraysArg)))
            fromDataArrays(:) = fromDataArraysArg(:)
          end if
        end if
        if (.not. allocated(fromDataArrays)) then
          log_debug(*) 'Data array "', trim(newTransferHelper%name), '" uses "fromFile" and not', &
                ' "fromDataArrays".'
          ! in order to have a properly parsed transfer function (in case it exists)
          ! the name of the TransferHelper is also the only data array used for the transfer function
          allocate(fromDataArrays(1))
          fromDataArrays(1) = newTransferHelper%name
          ! dummy allocation, so setting of transfer func works
          allocate(newTransferHelper%inputFieldNames(1))
          newTransferHelper%inputFieldNames(1) = newTransferHelper%name
        end if
      end if
    end if

    if (.not. newTransferHelper%doInitFromFile) then
      ! init from other DataArrays
      if (present(fromDataArraysArg)) then
        call newTransferHelper%check_and_set_input_field_names(fromDataArraysArg)
        allocate(fromDataArrays(size(fromDataArraysArg)))
        fromDataArrays(:) = fromDataArraysArg(:)
      else
        log_error(*) 'Data array "', trim(newTransferHelper%name), '" must be initialized either "fromFile" or', &
                ' "fromDataArrays", set this in the file', filenameNamelistMprDefault
        stop 1
      end if
    end if
    if (present(transferFuncName)) then
      if (trim(transferFuncName) /= trim(defaultAlias)) then
        if (present(fromParameters)) then
          if (size(fromParameters) > 0) then
            allocate(newTransferHelper%globalParameters(size(fromParameters)))
            newTransferHelper%globalParameters = fromParameters(:)
          end if
        end if
        ! check whether name exists
        call newTransferHelper%func_str_to_funcName(transferFuncName, fromDataArrays)
      else
        newTransferHelper%transferFuncName = defaultAlias
        allocate(newTransferHelper%globalParameters(0))
      end if
    else
      allocate(newTransferHelper%globalParameters(0))
      newTransferHelper%transferFuncName = defaultAlias
    end if
    log_subtrace(*) 'Set property transferFuncName to ', trim(newTransferHelper%transferFuncName)
    if (present(transferFuncLabel)) then
      newTransferHelper%transferFuncLabel = transferFuncLabel
    else
      newTransferHelper%transferFuncLabel = defaultAlias
    end if

    deallocate(fromDataArrays)
    log_subtrace(*) 'Finished Initialization of type newTransferHelper'

  end function newTransferHelper

  subroutine reset_transferHelper(self)
    class(TransferHelper), intent(inout) :: self

    self%name = ''
    self%id=0

    ! upscaling
    if (allocated(self%targetCoords)) deallocate(self%targetCoords)
    ! main
    if (allocated(self%inputFieldNames)) deallocate(self%inputFieldNames)
    if (allocated(self%globalParameters)) deallocate(self%globalParameters)

  end subroutine reset_transferHelper

  subroutine init_and_execute_transfer(self, data, reshapedMask, coords)
    class(TransferHelper), intent(inout) :: self
    real(dp), dimension(:), allocatable, intent(out) :: data
    logical, dimension(:), allocatable, intent(out) :: reshapedMask
    type(CoordinatePointer), allocatable, dimension(:), intent(out) :: coords

    if (self%doInitFromFile) then
      call self%read_nc(data, reshapedMask, coords)
    end if
    call self%execute_transfer(data, reshapedMask, coords)

  end subroutine init_and_execute_transfer

  subroutine check_and_set_input_field_names(self, inputFieldNames)
    class(TransferHelper), intent(inout) :: self
    character(*), intent(in), dimension(:) :: inputFieldNames

    integer(i4) :: i, n, iDA, nDA

    ! this allocates the vectors for the input field data and checks for correctness
    n = size(inputFieldNames)
    if (allocated(self%inputFieldNames)) deallocate(self%inputFieldNames)
    allocate(self%inputFieldNames(n))
    nDA = get_n_initialized(MPR_DATA_ARRAYS)

    do i = 1, n
      iDA = get_index_in_vector(inputFieldNames(i), MPR_DATA_ARRAYS)
      if (iDA > nDA .and. .not. (&
              index(trim(inputFieldNames(i)), upperBoundName) > 0 .or. &
              index(trim(inputFieldNames(i)), lowerBoundName) > 0) .or. &
              trim(inputFieldNames(i)) == trim(self%name)) then
        log_error(*)  "Data array '", trim(inputFieldNames(i)), "' is not found in available data arrays.", &
                " Add the data array in the file ", trim(filenameNamelistMprDefault), " for DataArray '", &
                trim(self%name), "'"
        stop 1
      end if
    end do
    self%inputFieldNames(:) = inputFieldNames(:)

  end subroutine check_and_set_input_field_names

  subroutine set_input_field_containers(self, inputFieldContainers)
    class(TransferHelper), intent(inout) :: self
    type(InputFieldContainer), allocatable, dimension(:), intent(inout) :: inputFieldContainers

    integer(i4) :: i, iDA

    allocate(inputFieldContainers(size(self%inputFieldNames)))
    ! this sets the pointers for the InputFieldContainer
    do i = 1, size(self%inputFieldNames)
      ! this checks whether the input_field_name is a special name and refers to a coordinate

      if (index(trim(self%inputFieldNames(i)), upperBoundName) > 0 .or. &
              index(trim(self%inputFieldNames(i)), lowerBoundName) > 0) then
        if (i == 1) then
          ! the input data arrays provided contain a coordinate related information as the first data array
          ! this is not okay
          log_error(*)  "Data array '", trim(self%inputFieldNames(i)), "' is referring to a coordinate property", &
                  " and must not be the first item of the data array list."
          stop 1
        end if
        ! cycle to the next input data array
        cycle
      end if
      ! set the input field container pointers
      iDA = get_index_in_vector(self%inputFieldNames(i), MPR_DATA_ARRAYS)
      inputFieldContainers(i) = InputFieldContainer(MPR_DATA_ARRAYS(iDA)%data, &
              MPR_DATA_ARRAYS(iDA)%reshapedMask, &
              MPR_DATA_ARRAYS(iDA)%coords)
    end do

  end subroutine set_input_field_containers

  subroutine call_transfer_func(self, inputFieldContainers, data)
    class(TransferHelper), intent(inout) :: self
    type(InputFieldContainer), allocatable, dimension(:), intent(inout) :: inputFieldContainers
    real(dp), dimension(:), allocatable, intent(inout), target :: data
    character(22) :: funcName

    ! select the transfer function from all known
    ! check here, if correct no of input fields and parameter is set for respective transfer function
    funcName = 'transfer_function_' // compress(num2str(mprTransferFunctionTable%get_index(self%transferFuncName)))
    log_debug(*) 'call_transfer_func: looked for ', trim(self%transferFuncName), ' in lookup table and got ', &
            trim(funcName)
    select case(trim(funcName))
    case('transfer_function_2')
      call self%check_transfer_func_args(1_i4, 0_i4, size(inputFieldContainers))
      data = transfer_function_2(inputFieldContainers, self%globalParameters)
    case('transfer_function_1')
      ! no transfer functions was specified -> set to identical mapping
      call self%check_transfer_func_args(1_i4, 0_i4, size(inputFieldContainers))
      data = transfer_function_1(inputFieldContainers, self%globalParameters)
    case default
      log_error(*) "1.) Please run 'python -m src_python.pre_proc.update_tfs_in_fortran_source' " // &
              "to add new transfer functions. Add the '--help' argument for details."
      log_error(*) "2.) Compile again."
      stop 1
    end select

  end subroutine call_transfer_func

  subroutine check_transfer_func_args(self, nArgs, nParameters, nFields)
    class(TransferHelper), intent(inout) :: self
    integer(i4), intent(in) :: nArgs
    integer(i4), intent(in) :: nParameters
    integer(i4), intent(in) :: nFields

    if (nFields /= nArgs) then
      log_error("(1X,A,A,A,A,A,A,I0,A,I0)")  "the transfer function '", trim(self%transferFuncName), &
              "' of effective parameter '", trim(self%name), "' received ", nFields, &
              " input fields, but needed ", nArgs
      stop 1
    end if
    if (nParameters > 0) then
      if (size(self%globalParameters) /= nParameters) then
        log_error("(1X,A,A,A,A,A,I0,A,I0)")  "the transfer function '", trim(self%transferFuncName), &
                "' of effective parameter '", &
                trim(self%name), "' received ", size(self%globalParameters), &
                " parameters, but needed ", nParameters
        stop 1
      end if
    end if

  end subroutine check_transfer_func_args

  subroutine execute_transfer(self, data, reshapedMask, coords)
    class(TransferHelper), intent(inout) :: self
    real(dp), dimension(:), allocatable, intent(inout), target :: data
    logical, dimension(:), allocatable, intent(inout), target :: reshapedMask
    type(CoordinatePointer), allocatable, dimension(:), intent(inout), target :: coords

    type(InputFieldContainer), allocatable, dimension(:) :: inputFieldContainers
    integer(i4) :: i

    if (.not. (allocated(data) .and. allocated(reshapedMask) .and. allocated(coords))) then
      ! check if all dimensions of input_fields are equal
      call self%check_data_arrays(reshapedMask, coords, inputFieldContainers)
    else
      ! we assume that data, mask and coords are only already allocated here, if read directly from nc file
      ! so we only have one input field container, which we now create and let point to allocated data
      allocate(inputFieldContainers(1))
      call add_temporal_data(data)
      i = get_temporal_data_index()

      ! this temporal data construct is needed, as the InputFieldContainer only contains pointers
      inputFieldContainers(1) = InputFieldContainer(MPR_TEMPORAL_DATA(i)%data, reshapedMask, coords)
    end if
    ! set and call the transfer_func (only here, because check_data_arrays optionally alters input_fields)
    ! this was previously done in here with a procedure pointer but was rebuilt due to compiler issues with Intel
    call self%call_transfer_func(inputFieldContainers, data)
    ! check for newly introduced missing values and apply to mask
    if (CHECK_FOR_NODATAVALUE) call adapt_mask(data, reshapedMask, coords)
    ! this clears temporal data (e.g. remasked input fields to make them inter-comparable) used for transfer function
    call clear_temporal_data()
    ! TODO: think about when it is safe to call clear_temporal_mask()
    deallocate(inputFieldContainers)

  end subroutine execute_transfer

  subroutine adapt_mask(data, reshapedMask, coords)
    real(dp), dimension(:), allocatable, intent(inout), target :: data
    logical, dimension(:), allocatable, intent(inout), target :: reshapedMask
    type(CoordinatePointer), allocatable, dimension(:), intent(in), target :: coords

    logical, dimension(:), allocatable :: temp_mask_1d
    logical, dimension(:, :), allocatable :: temp_mask_2d
    logical, dimension(:, :, :), allocatable :: temp_mask_3d
    logical, dimension(:, :, :, :), allocatable :: temp_mask_4d
    logical, dimension(:, :, :, :, :), allocatable :: temp_mask_5d
    real(dp), dimension(:), allocatable :: data_1d
    real(dp), dimension(:, :), allocatable :: data_2d
    real(dp), dimension(:, :, :), allocatable :: data_3d
    real(dp), dimension(:, :, :, :), allocatable :: data_4d
    real(dp), dimension(:, :, :, :, :), allocatable :: data_5d

    integer(i8), dimension(size(coords)) :: newShape
    integer(i4) :: iCoord

    ! check if there is a need to alter the mask
    ! TODO: this check is expensive, make it optional (maybe by global property)
    if (any(eq(data, nodata_dp))) then
      ! issue a warning
      log_warn(*) 'There are new missing values in the data array'
      ! get the correct shape of the data
      newShape = [(coords(iCoord)%coord_p%count, iCoord=1, size(coords))]
      select case(size(newShape(:)))
      case(1)
        temp_mask_1d = reshape(reshapedMask, newShape(1 : 1))
        ! unpack the data to the full shape
        data_1d = unpack(data, temp_mask_1d, nodata_dp)
        ! get the new mask
        temp_mask_1d = ne(data_1d, nodata_dp)
        ! pack it
        reshapedMask = pack(temp_mask_1d, .true.)
        ! pack the masked data
        data = pack(data_1d, temp_mask_1d)
        deallocate(temp_mask_1d, data_1d)
      case(2)
        temp_mask_2d = reshape(reshapedMask, newShape(1 : 2))
        ! unpack the data to the full shape
        data_2d = unpack(data, temp_mask_2d, nodata_dp)
        ! get the new mask
        temp_mask_2d = ne(data_2d, nodata_dp)
        ! pack it
        reshapedMask = pack(temp_mask_2d, .true.)
        ! pack the masked data
        data = pack(data_2d, temp_mask_2d)
        deallocate(temp_mask_2d, data_2d)
      case(3)
        temp_mask_3d = reshape(reshapedMask, newShape(1 : 3))
        ! unpack the data to the full shape
        data_3d = unpack(data, temp_mask_3d, nodata_dp)
        ! get the new mask
        temp_mask_3d = ne(data_3d, nodata_dp)
        ! pack it
        reshapedMask = pack(temp_mask_3d, .true.)
        ! pack the masked data
        data = pack(data_3d, temp_mask_3d)
        deallocate(temp_mask_3d, data_3d)
      case(4)
        temp_mask_4d = reshape(reshapedMask, newShape(1 : 4))
        ! unpack the data to the full shape
        data_4d = unpack(data, temp_mask_4d, nodata_dp)
        ! get the new mask
        temp_mask_4d = ne(data_4d, nodata_dp)
        ! pack it
        reshapedMask = pack(temp_mask_4d, .true.)
        ! pack the masked data
        data = pack(data_4d, temp_mask_4d)
        deallocate(temp_mask_4d, data_4d)
      case(5)
        temp_mask_5d = reshape(reshapedMask, newShape(1 : 5))
        ! unpack the data to the full shape
        data_5d = unpack(data, temp_mask_5d, nodata_dp)
        ! get the new mask
        temp_mask_5d = ne(data_5d, nodata_dp)
        ! pack it
        reshapedMask = pack(temp_mask_5d, .true.)
        ! pack the masked data
        data = pack(data_5d, temp_mask_5d)
        deallocate(temp_mask_5d, data_5d)
      end select
    end if

  end subroutine adapt_mask

  subroutine check_data_arrays(self, reshapedMask, coords, inputFieldContainers)
    class(TransferHelper), intent(inout) :: self
    logical, dimension(:), allocatable, intent(inout) :: reshapedMask
    type(CoordinatePointer), allocatable, dimension(:), intent(inout) :: coords
    type(InputFieldContainer), allocatable, dimension(:), intent(inout) :: inputFieldContainers

    type(CoordinatePointer), allocatable, dimension(:) :: IF1, IF2
    integer(i4) :: iIF, jIF, nIF, iCoord, nCoord, i, j
    logical :: startCorrect, endCorrect
    logical, dimension(:), allocatable :: packedMaskIIF, packedMaskJIF, tempMask
    real(dp), dimension(:), allocatable :: temporalData, temporalData1D
    real(dp), dimension(:,:), allocatable :: temporalData2D
    character(maxNameLength) :: errorMessageStem

    integer(i4) :: lowIndex, upIndex
    integer(i8) :: nTrailingCoord, nLeadingCoord
    integer(i4), dimension(:), allocatable :: coordIndex
    character(maxNameLength) :: dimName
    real(dp), dimension(:), allocatable :: values

    ! here the correct input data are pointed to
    call self%set_input_field_containers(inputFieldContainers)

    ! if concatenation occurs, the number of inputFieldContainers is reduced by 1
    ! thus we have to iteratively revisit this part
    iIF = 1
    concat_loop : do while (.true.)
      nIF = size(inputFieldContainers)
      ! get the number of Coordinates in the first referenced DataArray
      nCoord = size(inputFieldContainers(1)%coords_p)
      ! loop over possible input_field pairs, compare each with field 1
      IF1 = inputFieldContainers(iIF)%coords_p
      do jIF = iIF + 1, nIF
        ! if this field is referring to a coordinate property, cycle
        if (.not. associated(inputFieldContainers(jIF)%data_p)) then
          cycle
        end if
        IF2 = inputFieldContainers(jIF)%coords_p

        if (size(IF2) /= nCoord) then
          log_error("(1X,A,A,A,A,I0,A,A,A,I0,A)")  "the number of coordinates for input_field " , &
                  trim(self%inputFieldNames(jIF)) , " (" , size(IF2) , ") and input_field " , &
                  trim(self%inputFieldNames(iIF)) , " (" , nCoord , ") do not match. "
          stop 1
        end if

        log_debug(*) 'Comparing the coordinates of each predictor data array:'
        do iCoord = 1, size(IF1)
          log_debug(*) 'predictor: ', trim(self%inputFieldNames(iIF)), ', coordinate: ',trim(IF1(iCoord)%coord_p%name)
        end do
        do iCoord = 1, size(IF2)
          log_debug(*) 'predictor: ', trim(self%inputFieldNames(jIF)), ', coordinate: ', trim(IF2(iCoord)%coord_p%name)
        end do

        ! this block is checking each coordinate for its similarity and they are optionally concatenated
        do iCoord = 1, nCoord
          errorMessageStem = "the " // trim(num2str(iCoord)) // "th coordinate for input_field " // &
                  trim(self%inputFieldNames(iIF)) // ": '" // trim(IF1(iCoord)%coord_p%name) // "' and input_field " // &
                  trim(self%inputFieldNames(jIF)) // ": '" // trim(IF2(iCoord)%coord_p%name) // "' "
          ! check start and step
          startCorrect = eq(IF1(iCoord)%coord_p%bounds(1), IF2(iCoord)%coord_p%bounds(1))
          endCorrect = eq(IF1(iCoord)%coord_p%bounds(2), IF2(iCoord)%coord_p%bounds(2))

          if (startCorrect .and. endCorrect) then
            if (ne(IF1(iCoord)%coord_p%step, IF2(iCoord)%coord_p%step)) then
              call IF1(iCoord)%coord_p%get_stats()
              call IF2(iCoord)%coord_p%get_stats()
              log_error(*) trim(errorMessageStem), "have different step sizes."
              stop 1
            else if (IF1(iCoord)%coord_p%count /= IF2(iCoord)%coord_p%count) then
              call IF1(iCoord)%coord_p%get_stats()
              call IF2(iCoord)%coord_p%get_stats()
              log_error(*) trim(errorMessageStem), "have different number of entries."
              stop 1
            else if (.not. inputFieldContainers(1)%coords_p(iCoord)%coord_p == &
                    inputFieldContainers(iIF)%coords_p(iCoord)%coord_p) then
              call inputFieldContainers(1)%coords_p(iCoord)%coord_p%get_stats()
              call inputFieldContainers(iIF)%coords_p(iCoord)%coord_p%get_stats()
              log_error(*) trim(errorMessageStem), "are not equal."
              stop 1
            end if

          else if ((.not. startCorrect) .and. (.not. endCorrect)) then
            ! check whether start of the one is the end of the other
            if (eq(IF1(iCoord)%coord_p%bounds(1), IF2(iCoord)%coord_p%bounds(2))) then
              call self%concatenate_input_fields(inputFieldContainers, jIF, iIF, iCoord)
            else if (eq(IF1(iCoord)%coord_p%bounds(2), IF2(iCoord)%coord_p%bounds(1))) then
              call self%concatenate_input_fields(inputFieldContainers, iIF, jIF, iCoord)
            else
              call IF1(iCoord)%coord_p%get_stats()
              call IF2(iCoord)%coord_p%get_stats()
              log_error(*) trim(errorMessageStem), "are not equal and cannot be concatenated."
              stop 1
            end if
            ! the input data are modified, so we have to reenter the loop, as the jIF, nIF are not valid anymore
            cycle concat_loop
          else
            call IF1(iCoord)%coord_p%get_stats()
            call IF2(iCoord)%coord_p%get_stats()
            log_error(*) trim(errorMessageStem), "are not equal."
            stop 1
          end if
        end do

        ! this block is checking the masks for congruency

        ! the size of mask should be equal as the coords are already equal...
        ! so we now check the mask values for cases where they are not equal
        if (any(inputFieldContainers(iIF)%reshaped_mask_p .neqv. &
                inputFieldContainers(jIF)%reshaped_mask_p)) then

          ! now check for each others' congruency
          packedMaskIIF = pack(inputFieldContainers(jIF)%reshaped_mask_p, &
                  inputFieldContainers(iIF)%reshaped_mask_p)
          packedMaskJIF = pack(inputFieldContainers(iIF)%reshaped_mask_p, &
                  inputFieldContainers(jIF)%reshaped_mask_p)
          if (all(packedMaskJIF)) then
            ! this means that for all mask_jIF, there is also a valid iIF, thus we crop iIF to jIF
            temporalData = pack(inputFieldContainers(iIF)%data_p, packedMaskIIF)
            call add_temporal_data(temporalData)
            i = get_temporal_data_index()
            ! this temporal data construct is needed, as the InputFieldContainer only contains pointers
            inputFieldContainers(iIF) = InputFieldContainer(MPR_TEMPORAL_DATA(i)%data, &
                    inputFieldContainers(jIF)%reshaped_mask_p, &
                    inputFieldContainers(iIF)%coords_p)
            ! the data with iIF index (1) are modified, this is okay, no reentering the loop
            deallocate(packedMaskIIF, packedMaskJIF)
          else if (all(packedMaskIIF)) then
            ! this means that for all mask_iIF, there is also a valid jIF, thus we crop jIF to iIF
            temporalData = pack(inputFieldContainers(jIF)%data_p, packedMaskJIF)
            call add_temporal_data(temporalData)
            i = get_temporal_data_index()
            ! this temporal data construct is needed, as the InputFieldContainer only contains pointers
            inputFieldContainers(jIF) = InputFieldContainer(MPR_TEMPORAL_DATA(i)%data, &
                    inputFieldContainers(iIF)%reshaped_mask_p, &
                    inputFieldContainers(jIF)%coords_p)
            ! the input data are modified, so we have to reenter the loop, as the masks now have to be compared again
            deallocate(packedMaskIIF, packedMaskJIF)
            cycle concat_loop
          else
            log_warn(*)  "the masks for input_field " // &
                    trim(self%inputFieldNames(jIF)) // " and input_field " // &
              trim(self%inputFieldNames(iIF)) // " are not congruent. They are now cropped to shared mask."

            ! this means that for all mask_iIF, there is also a valid jIF, thus we crop jIF to iIF
            temporalData = pack(inputFieldContainers(jIF)%data_p, packedMaskJIF)
            call add_temporal_data(temporalData)
            i = get_temporal_data_index()

            tempMask = inputFieldContainers(iIF)%reshaped_mask_p .and. &
                inputFieldContainers(jIF)%reshaped_mask_p
            call add_temporal_mask(tempMask)
            j = get_temporal_mask_index()
            ! this temporal data construct is needed, as the InputFieldContainer only contains pointers
            inputFieldContainers(jIF) = InputFieldContainer(MPR_TEMPORAL_DATA(i)%data, &
                    MPR_TEMPORAL_MASK(j)%mask, &
                    inputFieldContainers(jIF)%coords_p)

            ! this means that for all mask_jIF, there is also a valid iIF, thus we crop iIF to jIF
            temporalData = pack(inputFieldContainers(iIF)%data_p, packedMaskIIF)
            call add_temporal_data(temporalData)
            i = get_temporal_data_index()
            ! this temporal data construct is needed, as the InputFieldContainer only contains pointers
            inputFieldContainers(iIF) = InputFieldContainer(MPR_TEMPORAL_DATA(i)%data, &
                    MPR_TEMPORAL_MASK(j)%mask, &
                    inputFieldContainers(iIF)%coords_p)
            deallocate(packedMaskIIF, packedMaskJIF)
            ! the input data are modified, so we have to reenter the loop, as the masks now have to be compared again
            cycle concat_loop
          end if
        end if

      end do
      ! successfully looped over all input fields and now we leave the outer loop
      exit concat_loop

    end do concat_loop

    ! all Coordinates and the mask are equal, use coords and mask of the 1st data array
    allocate(coords(nCoord))
    coords(:) = inputFieldContainers(1)%coords_p(:)
    reshapedMask = inputFieldContainers(1)%reshaped_mask_p

    ! now that everything is checked and set, init the dimension property input_fields
    do jIF = iIF + 1, nIF
      if (.not. associated(inputFieldContainers(jIF)%data_p)) then
        ! check which property shall be used
        dimName = self%inputFieldNames(jIF)
        upIndex = index(trim(dimName), upperBoundName)
        lowIndex = index(trim(dimName), lowerBoundName)
        if (upIndex > 0) then
          ! use the upper bound property
          ! get the dimension name whose property shall be used
          dimName = dimName(1 : upIndex - 1)
        else if (lowIndex > 0) then
          ! use the lower bound property
          ! get the dimension name whose property shall be used
          dimName = dimName(1 : lowIndex - 1)
        else
          log_error(*) 'Unexpected configuration, contact developer! (mo_mpr_data_array)'
          stop 1
        end if
        ! look up the coordinate name in coords also considering aliases
        coordIndex = get_index_in_coordinate(dimName, coords, checkAliasesArg = .true.)
        if (size(coordIndex) > 1_i4) then
          log_error(*) 'Unexpected configuration, contact developer! (mo_mpr_data_array). ', &
            'Concatenating coordinates with multiple elementary coordinate is not supported'
          stop 1
        end if
        if (coordIndex(1) > size(coords)) then
          log_error(*) 'The coordinate whose property should be used ("', trim(dimName), &
                  '") is not a valid coordinate for the input fields of data array: ', &
                  trim(self%name)
          stop 1
        end if
        allocate(values(coords(coordIndex(1))%coord_p%count))
        if (upIndex > 0) then
          values = coords(coordIndex(1))%coord_p%values(1 : coords(coordIndex(1))%coord_p%count)
        else if (lowIndex > 0) then
          values = coords(coordIndex(1))%coord_p%values(0 : coords(coordIndex(1))%coord_p%count - 1_i8)
        end if

        ! count the copies we have to pass on to the spread(values) command later on
        ! this depends on whether the coordinate is before coordIndex or after
        nTrailingCoord = 1_i8; nLeadingCoord = 1_i8
        do i = 1, size(coords)
          if (i < coordIndex(1)) then
            nLeadingCoord = nLeadingCoord * coords(i)%coord_p%count
          else if (i > coordIndex(1)) then
            nTrailingCoord = nTrailingCoord * coords(i)%coord_p%count
          end if
          log_debug(*) 'check_data_arrays: coordinate', i, nLeadingCoord, nTrailingCoord
        end do
        ! bring the values to the required shape and pack again
        ! use spread the values accordingly (only separate between coords before and after coordIndex)
        ! finally use mask to crop data
        ! TODO: how to make that quicker?
        temporalData2D = spread(values, 2, nTrailingCoord)
        temporalData1D = pack(temporalData2D, .true.)
        temporalData2D = spread(temporalData1D, 1, nLeadingCoord)
        temporalData1D = pack(temporalData2D, .true.)
        temporalData = pack(temporalData1D, reshapedMask)

        call add_temporal_data(temporalData)
        i = get_temporal_data_index()
        ! this temporal data construct is needed, as the InputFieldContainer only contains pointers
        inputFieldContainers(jIF) = InputFieldContainer(MPR_TEMPORAL_DATA(i)%data, &
                reshapedMask, &
                coords)
        deallocate(values, temporalData2D, temporalData1D)

      end if
    end do

  end subroutine check_data_arrays

  subroutine func_str_to_funcName(self, funcString, varNames)
    !< converts string to function name

    !< it deletes blanks and line breaks
    !< it replaces fromDataArrays by strings x1,...,xn
    !< it replaces fromParameters by strings p1,...,pn
    !< it replaces double blanks

    !< author: Stephan Thober
    !< created: Jun 2018

    ! input variables
    class(TransferHelper), intent(inout) :: self
    character(len = *), intent(in) :: funcString
    character(len = *), dimension(:), intent(in) :: varNames

    ! local variables
    character(len = maxStringLength) :: replaceString
    character(len = maxStringLength) :: deleteString
    integer(i4) :: ii, istr, stringLength, counter, currentIndex
    character(len = 5), dimension(38, 2) :: operatorLookup
    real(dp), dimension(:), allocatable :: tempParam
    integer(i4), dimension(:), allocatable :: tempIndex, valueSortedIndex
    character(maxNameLength), dimension(:), allocatable :: tempName

    ! setup operator dictionary for translation
    operatorLookup(1, :) = (/'    +', '   pl'/)
    operatorLookup(2, :) = (/'    -', '   mi'/)
    operatorLookup(3, :) = (/'   **', '   po'/)
    operatorLookup(4, :) = (/'    /', '   di'/)
    operatorLookup(5, :) = (/'    *', '   ti'/)
    operatorLookup(6, :) = (/'    (', '   bs'/)
    operatorLookup(7, :) = (/'    )', '   be'/)
    operatorLookup(8, :) = (/'  exp', '   ex'/)
    ! this is base 10 logarithm
    operatorLookup(9, :) = (/'log10', '   l1'/)
    ! this is the natural logarithm
    operatorLookup(10, :) = (/'  log', '   l2'/)
    operatorLookup(11, :) = (/' then', '   th'/)
    operatorLookup(12, :) = (/' else', '   el'/)
    operatorLookup(13, :) = (/'   if', '   if'/)
    operatorLookup(14, :) = (/'   <=', '   le'/)
    operatorLookup(15, :) = (/'    <', '   lt'/)
    operatorLookup(16, :) = (/'   >=', '   ge'/)
    operatorLookup(17, :) = (/'    >', '   gt'/)
    operatorLookup(18, :) = (/'   ==', '   eq'/)
    operatorLookup(19, :) = (/'.and.', '   ad'/)
    operatorLookup(20, :) = (/' .or.', '   or'/)
    operatorLookup(21, :) = (/'  abs', '   ab'/)
    operatorLookup(22, :) = (/' acos', '   ac'/)
    operatorLookup(23, :) = (/'  max', '   ma'/)
    operatorLookup(24, :) = (/'  min', '   mi'/)
    operatorLookup(25, :) = (/' asin', '   as'/)
    operatorLookup(26, :) = (/' atan', '   at'/)
    operatorLookup(27, :) = (/'atan2', '   au'/)
    operatorLookup(28, :) = (/' cosh', '   ch'/)
    operatorLookup(29, :) = (/' sinh', '   sh'/)
    operatorLookup(30, :) = (/' sqrt', '   sq'/)
    operatorLookup(31, :) = (/'    ^', '   po'/)
    operatorLookup(32, :) = (/'  end', '   en'/)
    operatorLookup(33, :) = (/'where', '   wh'/)
    operatorLookup(34, :) = (/'.not.', '   no'/)
    operatorLookup(35, :) = (/' tanh', '   tx'/)
    operatorLookup(36, :) = (/'  sin', '   si'/)
    operatorLookup(37, :) = (/'  cos', '   co'/)
    operatorLookup(38, :) = (/'  tan', '   ta'/)

    ! initialize
    replaceString = funcString
    deleteString = funcString
    log_subtrace(*) 'func_str_to_funcName: Setting property transferFuncName with ', trim(funcString)

    ! replace variable names
    do ii = 1, size(varNames)
      log_subtrace(*) 'func_str_to_funcName: replacing variable name ', trim(varNames(ii)), ' by ', 'x' // compress(num2str(ii))
      replaceString = replace_word(replaceString, trim(varNames(ii)), 'x' // compress(num2str(ii)), .true.)
      deleteString = replace_word(deleteString, trim(varNames(ii)), '', .true.)
    end do

    log_trace(*) 'func_str_to_funcName: replaced variable names in replaceString: ', trim(replaceString)
    log_trace(*) 'func_str_to_funcName: deleted variable names in deleteString: ', trim(deleteString)

    ! replace parameter names (only in case, there no parameters yet)
    if (.not. allocated(self%globalParameters)) then
      allocate(tempParam(size(MPR_PARAMETERS%names)))
      allocate(tempIndex(size(MPR_PARAMETERS%names)))
      allocate(tempName(size(MPR_PARAMETERS%names)))

      counter = 0
      ! loop over all existing parameters and check if they are in the replaceString
      do ii = 1, size(MPR_PARAMETERS%names)
        currentIndex = index_word(replaceString, trim(MPR_PARAMETERS%names(ii)), .true.)
        if (currentIndex > 0) then
          counter = counter + 1
          ! add the parameter value to temporary vector of parameters
          tempParam(counter) = MPR_PARAMETERS%values(ii)
          tempIndex(counter) = currentIndex
          tempName(counter) = MPR_PARAMETERS%names(ii)
        end if
      end do
      ! sort all found items by index, so we always get the same order of p1, p2, p3 in replaceString no matter
      ! their order in the nml
      allocate(valueSortedIndex(counter))

      ! get index of sort values
      valueSortedIndex(:) = sort_index(tempIndex(1 : counter))

      allocate(self%globalParameters(counter))
      do ii = 1, counter
        replaceString = replace_word(replaceString, compress(tempName(valueSortedIndex(ii))), &
                'p' // compress(num2str(ii)), .true.)
        deleteString = replace_word(deleteString, compress(tempName(valueSortedIndex(ii))), '', .true.)
        self%globalParameters(ii) = tempParam(valueSortedIndex(ii))
      end do
      ! move the values from temporary to final vector
      deallocate(tempParam, tempIndex, tempName, valueSortedIndex)
    end if

    log_trace(*) 'func_str_to_funcName: replaced parameter names in replaceString: ', trim(replaceString)
    log_trace(*) 'func_str_to_funcName: deleted parameter names in deleteString: ', trim(deleteString)

    ! replace all operators
    do ii = 1, size(operatorLookup, dim = 1)
      replaceString = Replace_Text(replaceString, compress(operatorLookup(ii, 1)), '_' // trim(compress(operatorLookup(ii, 2))) // '_')
      deleteString = Replace_Text(deleteString, compress(operatorLookup(ii, 1)), '')
    end do

    log_trace(*) 'func_str_to_funcName: replaced operators in replaceString: ', trim(replaceString)
    log_trace(*) 'func_str_to_funcName: deleted operators in deleteString: ', trim(deleteString)

    ! delete parameters
    do ii = maxNoParameterPerDA, 1, -1
      deleteString = Replace_Text(deleteString, 'p' // compress(num2str(ii)), '')
    end do

    ! eliminate blanks
    replaceString = compress(replaceString)
    deleteString = compress(deleteString)

    ! remove possible underscores at start and end
    stringLength = len_trim(replaceString)
    istr = scan(replaceString, '_', .false.)
    if (istr == 1_i4) then
      replaceString = trim(replaceString(2 : stringLength))
      stringLength = len_trim(replaceString)
    end if
    istr = scan(replaceString, '_', .true.)
    if (istr == stringLength) then
      replaceString = trim(replaceString(1 : stringLength - 1))
    end if
    ! remove duplicates
    replaceString = replace_text(replaceString, '__', '_')

    if (trim(deleteString) /= '') then
      log_warn(*) "could not correctly parse the string for the transfer function: " // trim(funcString)
      log_warn(*) "  those chars could not be matched to DataArrays, Parameters or operators: " // trim(deleteString)
      replaceString = funcString
    else
      log_debug(*) 'translated string: ' // trim(funcString) // ' -> ' // trim(replaceString)
    end if

    self%transferFuncName = replaceString

  end subroutine func_str_to_funcName

  subroutine get_stats(self)
    !< prints some information on the TransferHelper
    !< useful for debugging or printing before raising error messages

    !> Coordinate type-bound procedure
    class(TransferHelper), intent(in) :: self

    integer(i4) :: iInputField

    log_info(*) 'TransferHelper "', trim(self%name), '" with id ', self%id
    if (self%doInitFromFile) then
      log_info(*) 'It is initialized from file "', trim(self%fileName), '" using transfer function:', &
              trim(self%transferFuncName)
    else
      log_info(*) 'It is initialized from other data arrays using transfer function:', trim(self%transferFuncName)
    end if

    log_info(*) 'It uses the data arrays:'
    do iInputField=1, size(self%inputFieldNames)
      log_info(*) iInputField, ': ', trim(self%inputFieldNames(iInputField))
    end do

    log_info(*) 'It uses the parameters:'
    do iInputField=1, size(self%globalParameters)
      log_info(*) iInputField, ':', self%globalParameters(iInputField)
    end do

  end subroutine get_stats


  subroutine read_nc(self, data, reshapedMask, coords)
    class(TransferHelper), intent(inout) :: self
    real(dp), dimension(:), allocatable, intent(out) :: data
    logical, dimension(:), allocatable, intent(out) :: reshapedMask
    type(CoordinatePointer), allocatable, dimension(:), intent(out) :: coords

    type(NcDataset) :: nc
    type(NcVariable) :: ncVar
    type(NcDimension), dimension(:), allocatable :: NcDims
    integer(i4), dimension(:), allocatable :: varShape

    character(maxNameLength), allocatable, dimension(:) :: coordNames
    integer(i4), dimension(1) :: id
    integer(i4) :: i, iCoord, nCoords
    character(6) :: stagger
    real(dp) :: bound
    type(Coordinate) :: coord_

    ! read Dataset and store the variable of interest
    log_debug(*) 'read_nc: reading DataArray from file ', trim(self%fileName)
    nc = NcDataset(self%fileName, "r")
    ncVar = nc%getVariable(trim(self%name))

    ! infer the dimensions...
    varShape = ncVar%getShape()
    ! then send these to the datagetter and packer
    select case(size(varShape))
    case(1)
      call get_mask_and_packed_data1D(ncVar, reshapedMask, data)
    case(2)
      call get_mask_and_packed_data2D(ncVar, reshapedMask, data)
    case(3)
      call get_mask_and_packed_data3D(ncVar, reshapedMask, data)
    case(4)
      call get_mask_and_packed_data4D(ncVar, reshapedMask, data)
    case default
      log_error("(1X,A,I0,A)") size(varShape), "-dimensional arrays are not supported as input fields."
      stop 1
    end select

    ! now get the dimensions
    allocate(coordNames(size(varShape)))
    allocate(coords(size(varShape)))
    ! get all dimension names of the variable
    NcDims = ncVar%getDimensions()
    do i = 1, size(varShape)
      coordNames(i) = trim(NcDims(i)%getName())
    end do

    ! close the file, it will be reopened during Coordinate initialization
    call nc%close()

    log_debug(*) 'read_nc: read DataArray "', trim(self%name), '" with shape ', varShape, ' and packed size ', &
          size(data, kind=i8)

    do iCoord = 1, size(varShape)
      ! set the stagger value
      stagger = trim(defaultCoordStagger)
      bound = nodata_dp

      ! try to override default if already in dimensions
      id = get_index_in_coordinate(coordNames(iCoord))
      nCoords = get_n_initialized(MPR_COORDINATES)
      if (id(1) <= nCoords) then
        ! use settings provided in namelist and update the dimension at index id(1)
        stagger = MPR_COORDINATES(id(1))%staggerName
        bound = MPR_COORDINATES(id(1))%get_bound()
      end if
      coord_ = Coordinate(&
              name=coordNames(iCoord), &
              id=id(1), &
              stagger=stagger, &
              fileName=self%fileName, &
              bound=bound &
      )
      ! only update, if not yet finalized
      if (id(1) > nCoords) then
        call add_coordinate(coord_, id(1))
      else if (.not. MPR_COORDINATES(id(1))%is_finalized()) then
        call add_coordinate(coord_, id(1))
        ! if coordinate is finalized, but the vector is different, add a new coordinate with same name
      elseif (coord_%is_finalized() .and. MPR_COORDINATES(id(1)) /= coord_) then
        id(1) = nCoords + 1
        coord_%id = id(1)
        call add_coordinate(coord_, id(1))
      end if
      ! and now set the pointer correctly
      call coords(iCoord)%set_coordinate_pointer(id(1))

    end do

  end subroutine read_nc

  subroutine concatenate_input_fields(self, IF, iIF, jIF, coordIndex)
    !> the TransferHelper class
    class(TransferHelper), intent(inout) :: self
    !> all the InputFieldContainers of the current TransferHelper
    type(InputFieldContainer), dimension(:), allocatable, intent(inout) :: IF
    !> index of first DataArray to be concatenated
    integer(i4), intent(in) :: iIF
    !> index of second DataArray to be concatenated
    integer(i4), intent(in) :: jIF
    !> index of Coordinate to concatenate along
    integer(i4), intent(in) :: coordIndex

    integer(i4) :: i, j, id
    character(maxNameLength), allocatable, dimension(:) :: newInputFieldNames
    character(maxNameLength) :: newDataArrayName
    real(dp), dimension(:), allocatable :: newData
    logical, dimension(:), allocatable :: newReshapedMask
    type(CoordinatePointer), dimension(:), allocatable :: newCoords
    type(DataArray) :: dataarray_

    newDataArrayName = trim(self%inputFieldNames(iIF)) // '_cat_' // trim(self%inputFieldNames(jIF))

    log_debug(*) ' For data array: ', trim(self%name)
    log_debug(*) ' concatenating data arrays'
    log_debug(*) ' field 1:   ', trim(self%inputFieldNames(iIF))
    log_debug(*) ' field 2:   ', trim(self%inputFieldNames(jIF))
    log_debug(*) ' to temporary new field: ', trim(newDataArrayName)

    ! create the new data, mask and coords (also new Coordinate) from two InputFieldContainers
    call IF(iIF)%concat(other = IF(jIF), &
            coordIndex = coordIndex, &
            data = newData, reshapedMask = newReshapedMask, coords = newCoords)

    ! --- create new data array
    ! first get the index at which it should be inserted in global vector
    id = get_index_in_vector(newDataArrayName, MPR_DATA_ARRAYS)
    ! ... init the derived type ...
    dataarray_ = DataArray(&
            name = newDataArrayName, &
            id = id, &
            fromData = newData, &
            fromReshapedMask = newReshapedMask, &
            fromCoords = newCoords, &
            toFileArg = .false. &
            )

    ! ... and add it to global vector
    call add_data_array(dataarray_, id)

    ! --- update self
    allocate(newInputFieldNames(size(self%inputFieldNames) - 1_i4))

    j = 1
    do i = 1, size(self%inputFieldNames)
      if (i == jIF) then
        ! skip if we are at index of dimension that is removed
        cycle
      end if
      if (i == iIF) then
        ! set new dimension name if we are at index of dimension that is updated
        newInputFieldNames(j) = newDataArrayName
      else
        newInputFieldNames(j) = self%inputFieldNames(i)
      end if
      j = j + 1
    end do

    ! set the new_input_data
    call self%check_and_set_input_field_names(newInputFieldNames(:))
    deallocate(IF)
    call self%set_input_field_containers(IF)

    ! --- clean up
    deallocate(newData, newReshapedMask, newCoords, newInputFieldNames)

  end subroutine concatenate_input_fields


end module mo_mpr_data_array

