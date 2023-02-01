#include "flogging.h"
module mo_mpr_input_field_container

  use mo_kind, only : dp, i4, i8
  use mo_mpr_coordinate, only: CoordinatePointer, Coordinate, MPR_COORDINATES, add_coordinate, get_common_unit
  use mo_constants, only: nodata_dp
  use mo_mpr_constants, only: maxNameLength, maxStringLength
  use mo_mpr_utils, only: get_index_in_vector, get_n_initialized
  use mo_append, only: append
  use mo_mpr_util_types, only: mpr_add_coord_alias
  use flogging

  implicit none

  private

  ! --------------------------------------------------------------------------------------
  type, public :: InputFieldContainer

    real(dp), dimension(:), pointer, public :: data_p => null()
    logical, dimension(:), pointer, public :: reshaped_mask_p => null()
    type(CoordinatePointer), dimension(:), pointer, public :: coords_p => null()

  contains
    private
    procedure, public :: concat
    procedure :: append_input_field_data
    procedure :: append_input_field_data_1d
    procedure :: append_input_field_data_2d
    procedure :: append_input_field_data_3d
    procedure :: append_input_field_data_4d


  end type InputFieldContainer

  interface InputFieldContainer
    procedure newInputFieldContainer
  end interface InputFieldContainer

contains

  function newInputFieldContainer(data, reshapedMask, coords)
    real(dp), dimension(:), intent(in), target :: data
    logical, dimension(:), intent(in), target :: reshapedMask
    type(CoordinatePointer), dimension(:), intent(in), target :: coords
    type(InputFieldContainer):: newInputFieldContainer

    newInputFieldContainer%data_p => data
    newInputFieldContainer%reshaped_mask_p => reshapedMask
    newInputFieldContainer%coords_p => coords

  end function newInputFieldContainer

  subroutine concat(self, other, coordIndex, data, reshapedMask, coords)
    class(InputFieldContainer), intent(inout):: self
    type(InputFieldContainer), intent(in):: other
    integer(i4), intent(in):: coordIndex
    real(dp), dimension(:), allocatable, intent(out) :: data
    logical, dimension(:), allocatable, intent(out) :: reshapedMask
    type(CoordinatePointer), dimension(:), allocatable, intent(out) :: coords

    integer(i8), allocatable :: dataShape(:, :)
    real(dp), allocatable :: coordValues(:)
    character(maxNameLength) :: newCoordinateName
    character(maxStringLength) :: unit
    integer(i4) :: iCoord, nCoord, id, nCoordsGlobal
    type(Coordinate) :: coord_

    ! -- create the new coordinate
    ! create the new name
    newCoordinateName = trim(self%coords_p(coordIndex)%coord_p%name) // &
            '_cat_' // trim(other%coords_p(coordIndex)%coord_p%name)
    ! create the new values
    coordValues = self%coords_p(coordIndex)%coord_p%values(:)
    call append(coordValues, other%coords_p(coordIndex)%coord_p%values(1:))
    unit = get_common_unit(self%coords_p(coordIndex)%coord_p, other%coords_p(coordIndex)%coord_p)

    ! first get the index at which it should be inserted in global vector
    nCoordsGlobal = get_n_initialized(MPR_COORDINATES)
    id = get_index_in_vector(newCoordinateName, MPR_COORDINATES)
    ! ... init the derived type ...
    log_trace(*) 'concat: initializing the coordinate ', trim(newCoordinateName)
    coord_ = Coordinate(&
            name=newCoordinateName, &
            id=id, &
            stagger=self%coords_p(coordIndex)%coord_p%staggerName, &
            values=coordValues(2:), &
            bound=self%coords_p(coordIndex)%coord_p%bounds(1), &
            unit=unit &
            )
    ! --- update the aliases
    call mpr_add_coord_alias(newCoordinateName, trim(self%coords_p(coordIndex)%coord_p%name))

    ! --- update MPR_COORDINATES
    ! only update, if not yet finalized
    if (.not. MPR_COORDINATES(id)%is_finalized()) then
      log_trace(*) 'concat: updating the coordinate ', trim(newCoordinateName), ' in global list'
      call add_coordinate(coord_, id)
    ! if dim is finalized, but the values is different, add a new coordinate with same name
    elseif (MPR_COORDINATES(id) /= coord_) then
      log_trace(*) 'concat: adding the coordinate ', trim(newCoordinateName), ' to global list'
      id = nCoordsGlobal + 1
      coord_%id = id
      call add_coordinate(coord_, id)
    end if

    nCoord = size(self%coords_p)
    allocate(coords(nCoord))
    allocate(dataShape(2, nCoord))

    do iCoord = 1, nCoord
      ! create the new Coordinate pointers
      if (iCoord == coordIndex) then
        call coords(iCoord)%set_coordinate_pointer(id)
      else
        call coords(iCoord)%set_coordinate_pointer(self%coords_p(iCoord)%coord_p%id)
      end if
      ! --- append input field data
      dataShape(1, iCoord) = self%coords_p(iCoord)%coord_p%count
      dataShape(2, iCoord) = other%coords_p(iCoord)%coord_p%count
    end do

    call self%append_input_field_data(dataShape(1, :), dataShape(2, :), other, coordIndex, &
            data, reshapedMask)

  end subroutine concat

  subroutine append_input_field_data(self, dataShape1, dataShape2, other, coordIndex, data, reshapedMask)
    class(InputFieldContainer), intent(in) :: self
    integer(i8) :: dataShape1(:), dataShape2(:)
    type(InputFieldContainer), intent(in) :: other
    integer(i4), intent(in) :: coordIndex
    real(dp), dimension(:), allocatable, intent(out) :: data
    logical, dimension(:), allocatable, intent(out) :: reshapedMask

    select case(size(dataShape1))
    case(1)
      call self%append_input_field_data_1d(other, coordIndex, data, reshapedMask)
    case(2)
      call self%append_input_field_data_2d(dataShape1, dataShape2, other, coordIndex, data, reshapedMask)
    case(3)
      call self%append_input_field_data_3d(dataShape1, dataShape2, other, coordIndex, data, reshapedMask)
    case(4)
      call self%append_input_field_data_4d(dataShape1, dataShape2, other, coordIndex, data, reshapedMask)
    end select

  end subroutine append_input_field_data

  subroutine append_input_field_data_1d(self, other, coordIndex, data, reshapedMask)

    class(InputFieldContainer), intent(in) :: self
    type(InputFieldContainer), intent(in) :: other
    integer(i4), intent(in) :: coordIndex
    real(dp), dimension(:), allocatable, intent(out) :: data
    logical, dimension(:), allocatable, intent(out) :: reshapedMask

    if (coordIndex /= 1) then
      log_error("(1X,A,A,I0,A)") 'append_input_field_data_1d: coordIndex is ', coordIndex, &
      ' but has to be 1 for 1-dimensional array '
      stop 1
    end if

    reshapedMask = self%reshaped_mask_p
    call append(reshapedMask, other%reshaped_mask_p)

    data = self%data_p
    call append(data, other%data_p)

  end subroutine append_input_field_data_1d

  subroutine append_input_field_data_2d(self, dataShape1, dataShape2, other, coordIndex, data, reshapedMask)

    class(InputFieldContainer), intent(in) :: self
    integer(i8) :: dataShape1(:), dataShape2(:)
    type(InputFieldContainer), intent(in) :: other
    integer(i4), intent(in) :: coordIndex
    real(dp), dimension(:), allocatable, intent(out) :: data
    logical, dimension(:), allocatable, intent(out) :: reshapedMask

    logical, allocatable :: m1(:, :), m2(:, :) ! temporary masks
    real(dp), allocatable :: d1(:, :), d2(:, :) ! temporary data

    m1 = reshape(self%reshaped_mask_p, dataShape1(1:2))
    m2 = reshape(other%reshaped_mask_p, dataShape2(1:2))

    d1 = unpack(self%data_p, mask=m1, field=nodata_dp)
    d2 = unpack(other%data_p, mask=m2, field=nodata_dp)

    call append(m1, m2, idim=coordIndex)
    call append(d1, d2, idim=coordIndex)

    reshapedMask = reshape(m1, (/product(shape(m1))/))
    data = pack(d1, mask=m1)

    deallocate(m1, m2, d1, d2)

  end subroutine append_input_field_data_2d

  subroutine append_input_field_data_3d(self, dataShape1, dataShape2, other, coordIndex, data, reshapedMask)

    class(InputFieldContainer), intent(in) :: self
    integer(i8) :: dataShape1(:), dataShape2(:)
    type(InputFieldContainer), intent(in) :: other
    integer(i4), intent(in) :: coordIndex
    real(dp), dimension(:), allocatable, intent(out) :: data
    logical, dimension(:), allocatable, intent(out) :: reshapedMask

    logical, allocatable :: m1(:, :, :), m2(:, :, :) ! temporary masks
    real(dp), allocatable :: d1(:, :, :), d2(:, :, :) ! temporary data

    m1 = reshape(self%reshaped_mask_p, dataShape1(1:3))
    m2 = reshape(other%reshaped_mask_p, dataShape2(1:3))

    d1 = unpack(self%data_p, mask=m1, field=nodata_dp)
    d2 = unpack(other%data_p, mask=m2, field=nodata_dp)

    call append(m1, m2, idim=coordIndex)
    call append(d1, d2, idim=coordIndex)

    reshapedMask = reshape(m1, (/product(shape(m1))/))
    data = pack(d1, mask=m1)

    deallocate(m1, m2, d1, d2)

  end subroutine append_input_field_data_3d

  subroutine append_input_field_data_4d(self, dataShape1, dataShape2, other, coordIndex, data, reshapedMask)

    class(InputFieldContainer), intent(in) :: self
    integer(i8) :: dataShape1(:), dataShape2(:)
    type(InputFieldContainer), intent(in) :: other
    integer(i4), intent(in) :: coordIndex
    real(dp), dimension(:), allocatable, intent(out) :: data
    logical, dimension(:), allocatable, intent(out) :: reshapedMask

    logical, allocatable :: m1(:, :, :, :), m2(:, :, :, :) ! temporary masks
    real(dp), allocatable :: d1(:, :, :, :), d2(:, :, :, :) ! temporary data

    m1 = reshape(self%reshaped_mask_p, dataShape1(1:4))
    m2 = reshape(other%reshaped_mask_p, dataShape2(1:4))

    d1 = unpack(self%data_p, mask=m1, field=nodata_dp)
    d2 = unpack(other%data_p, mask=m2, field=nodata_dp)

    call append(m1, m2, idim=coordIndex)
    call append(d1, d2, idim=coordIndex)

    reshapedMask = reshape(m1, (/product(shape(m1))/))
    data = pack(d1, mask=m1)

    deallocate(m1, m2, d1, d2)

  end subroutine append_input_field_data_4d


end module mo_mpr_input_field_container
