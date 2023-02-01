!> \file mo_mpr_utils.f90

!> \brief MPR related types

!> \details A wrapper around mo_netcdf
!
!> \authors Robert Schweppe
!> \date Jan 2018
#include "flogging.h"
module mo_mpr_utils

  use mo_kind, only: i4, dp, i8
  use mo_mpr_util_types, only: MprBaseType
  use flogging
  use mo_netcdf, only : NcDataset, NcDimension, NcVariable

  implicit none

  private

  public :: get_index_in_vector, get_n_initialized, read_coordinate_from_file, read_variable_from_file, &
          get_mask_and_packed_data1D, &
          get_mask_and_packed_data2D, &
          get_mask_and_packed_data3D, &
          get_mask_and_packed_data4D

  interface read_variable_from_file
    procedure read_variable_from_file_2d_dp, read_variable_from_file_1d_dp, &
            read_variable_from_file_1d_i4
  end interface read_variable_from_file

  ! --------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------
contains

  function get_index_in_vector(itemName, globalVector) result(index)
    !< returns the index of the MprBaseType of globalVector if it matches itemName, else the first free index

    !> the name of the MprBaseType to be found
    character(*), intent(in) :: itemName
    !> the array of MprBaseTypes to be searched
    class(MprBaseType), intent(in), dimension(:) :: globalVector
    integer(i4) :: index

    integer(i4) :: nItems

    ! get length of MPR_COORDINATES
    nItems = get_n_initialized(globalVector)
    ! check if item is already contained in coordinates
    do index = 1, nItems
      if (trim(globalVector(index)%name) == trim(itemName)) then
        ! if already contained, log information and return
        log_trace(*) 'get_index_in_vector: Item ',trim(itemName), ' already contained in search list'
        return
      end if
    end do
    log_trace(*) 'get_index_in_vector: Item ',trim(itemName), ' not contained in search list'

  end function get_index_in_vector

  function get_n_initialized(globalVector) result(n)
    !> returns the number of the MprBaseType items in globalVector that are initialized
    class(MprBaseType), intent(in), dimension(:) :: globalVector
    integer(i4) :: n
    integer(i4) :: iItem
    ! check for proper name
    n = 0_i4
    do iItem = 1, size(globalVector)
      if (globalVector(iItem)%is_initialized) then
        n = iItem
      end if
    end do

  end function get_n_initialized

  subroutine read_coordinate_from_file(nc, name, property)
    type(NcDataset), intent(inout) :: nc
    character(*), intent(in) :: name
    integer(i4), intent(out) :: property

    type(NcDimension) :: NcDim

    if (nc%hasDimension(trim(name))) then
      NcDim = nc%getDimension(trim(name))
      property = NcDim%getLength()
    else
      log_error(*) 'The coordinate "', trim(name), '" in source file ', trim(nc%fname), ' cannot be read'
      call nc%close()
      stop 1
    end if

  end subroutine read_coordinate_from_file

  subroutine read_variable_from_file_2d_dp(nc, name, property, mask)
    type(NcDataset), intent(inout) :: nc
    character(*), intent(in) :: name
    real(dp), dimension(:,:), allocatable, intent(out) :: property
    logical, dimension(:,:), allocatable, intent(out), optional :: mask

    type(NcVariable) :: ncVar
    real(dp) :: noDataValue, minValue, maxValue, scale, offset
    logical :: status

    if (nc%hasVariable(trim(name))) then
      ncVar = nc%getVariable(trim(name))
      call ncVar%getData(property, mask=mask)
    else
      log_error(*) 'The variable "', trim(name), '" in source file ', trim(nc%fname), ' cannot be read'
      call nc%close()
      stop 1
    end if
    call get_scale_and_offset(ncVar, scale, offset, status)
    if (status) then
      property = property * scale + offset
    end if

  end subroutine read_variable_from_file_2d_dp

  subroutine read_variable_from_file_1d_dp(nc, name, property, mask)
    type(NcDataset), intent(inout) :: nc
    character(*), intent(in) :: name
    real(dp), dimension(:), allocatable, intent(out) :: property
    logical, dimension(:), allocatable, intent(out), optional :: mask

    type(NcVariable) :: ncVar
    real(dp) :: noDataValue, minValue, maxValue, scale, offset
    logical :: status

    if (nc%hasVariable(trim(name))) then
      ncVar = nc%getVariable(trim(name))
      call ncVar%getData(property, mask=mask)
    else
      log_error(*) 'The variable "', trim(name), '" in source file ', trim(nc%fname), ' cannot be read'
      call nc%close()
      stop 1
    end if
    call get_scale_and_offset(ncVar, scale, offset, status)
    if (status) then
      property = property * scale + offset
    end if

  end subroutine read_variable_from_file_1d_dp

  subroutine read_variable_from_file_1d_i4(nc, name, property, mask)
    type(NcDataset), intent(inout) :: nc
    character(*), intent(in) :: name
    integer(i4), dimension(:), allocatable, intent(out) :: property
    logical, dimension(:), allocatable, intent(out), optional :: mask

    type(NcVariable) :: ncVar

    if (nc%hasVariable(trim(name))) then
      ncVar = nc%getVariable(trim(name))
      call ncVar%getData(property, mask=mask)
    else
      log_error(*) 'The variable "', trim(name), '" in source file ', trim(nc%fname), ' cannot be read'
      call nc%close()
      stop 1
    end if

  end subroutine read_variable_from_file_1d_i4

  subroutine get_mask_and_packed_data1D(ncVar, packMask, packedData)
    class(NcVariable), intent(inout) :: ncVar
    logical, dimension(:), allocatable, intent(out) :: packMask
    real(dp), dimension(:), allocatable, intent(out) :: packedData

    real(dp), dimension(:), allocatable :: ncData
    logical, dimension(:), allocatable :: ncMask
    real(dp) :: noDataValue, minValue, maxValue, scale, offset
    logical :: status

    ! get data
    call ncVar%getData(ncData, mask=ncMask)
    packMask = pack(ncMask, .true.)
    packedData = pack(ncData, ncMask)
    deallocate(ncData)
    deallocate(ncMask)
    call get_scale_and_offset(ncVar, scale, offset, status)
    if (status) then
      packedData = packedData * scale + offset
    end if

  end subroutine get_mask_and_packed_data1D

  subroutine get_mask_and_packed_data2D(ncVar, packMask, packedData)
    class(NcVariable), intent(inout) :: ncVar
    logical, dimension(:), allocatable, intent(out) :: packMask
    real(dp), dimension(:), allocatable, intent(out) :: packedData
    real(dp), dimension(:, :), allocatable :: ncData
    logical, dimension(:, :), allocatable :: ncMask
    real(dp) :: noDataValue, minValue, maxValue, scale, offset
    logical :: status

    ! get data
    call ncVar%getData(ncData, mask=ncMask)
    packMask = pack(ncMask, .true.)
    packedData = pack(ncData, ncMask)
    deallocate(ncData)
    deallocate(ncMask)
    call get_scale_and_offset(ncVar, scale, offset, status)
    if (status) then
      packedData = packedData * scale + offset
    end if

  end subroutine get_mask_and_packed_data2D

  subroutine get_mask_and_packed_data3D(ncVar, packMask, packedData)
    class(NcVariable), intent(inout) :: ncVar
    logical, dimension(:), allocatable, intent(out) :: packMask
    real(dp), dimension(:), allocatable, intent(out) :: packedData
    real(dp), dimension(:, :, :), allocatable :: ncData
    logical, dimension(:, :, :), allocatable :: ncMask
    real(dp) :: noDataValue, minValue, maxValue, scale, offset
    logical :: status

    ! get data
    call ncVar%getData(ncData, mask=ncMask)
    packMask = pack(ncMask, .true.)
    packedData = pack(ncData, ncMask)
    deallocate(ncData)
    deallocate(ncMask)
    call get_scale_and_offset(ncVar, scale, offset, status)
    if (status) then
      packedData = packedData * scale + offset
    end if


  end subroutine get_mask_and_packed_data3D

  subroutine get_mask_and_packed_data4D(ncVar, packMask, packedData)
    class(NcVariable), intent(inout) :: ncVar
    logical, dimension(:), allocatable, intent(out) :: packMask
    real(dp), dimension(:), allocatable, intent(out) :: packedData
    real(dp), dimension(:, :, :, :), allocatable :: ncData
    logical, dimension(:, :, :, :), allocatable :: ncMask

    real(dp) :: noDataValue, minValue, maxValue, scale, offset
    logical :: status

    ! get data
    call ncVar%getData(ncData, mask=ncMask)
    packMask = pack(ncMask, .true.)
    packedData = pack(ncData, ncMask)
    deallocate(ncData)
    deallocate(ncMask)
    call get_scale_and_offset(ncVar, scale, offset, status)
    if (status) then
      packedData = packedData * scale + offset
    end if

  end subroutine get_mask_and_packed_data4D

  subroutine get_scale_and_offset(ncVar, scale, offset, status)
    !< check whether scale and add_offset is set for this variable
    !> ncVariable to be checked for attributes
    type(NcVariable), intent(inout) :: ncVar
    !> value indicating a factor to be applied
    real(dp), intent(out) :: scale
    !> value indicating a summand to be applied
    real(dp), intent(out) :: offset
    !> boolean indicating whether values were found
    logical, intent(out) :: status

    ! determine values
    status = .false.
    scale = 1.0_dp
    offset = 0.0_dp
    if (ncVar%hasAttribute("scale_factor")) then
      call ncVar%getAttribute('scale_factor', scale)
      status = .true.
    end if

    if (ncVar%hasAttribute("add_offset")) then
      call ncVar%getAttribute('add_offset', offset)
      status = .true.
    end if

  end subroutine get_scale_and_offset

end module mo_mpr_utils

