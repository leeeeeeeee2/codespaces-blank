!> \file mo_mpr_coordinate.f90

!> \brief MPR related types

!> \details Specifications of MPR coordinates
!
!> \authors Robert Schweppe, Stephan Thober
!> \date Jan 2018

#include "flogging.h"
module mo_mpr_coordinate

  use mo_kind, only : sp, dp, i1, i2, i4, i8
  use mo_string_utils, only : num2str, compress
  use mo_constants, only : nodata_dp, nodata_i4, nodata_i8
  use mo_mpr_constants, only : defaultAlias, MaxNoCoords, maxTolerance, maxNameLength, maxStringLength, &
          defaultCoordUnits, defaultCoordIsAscending, defaultCoordCount, maxNoAttributes, maxNoAttributeValues
  use mo_utils, only : eq, ne
  use mo_netcdf, only : NcDataset, NcDimension, NcVariable
  use ieee_arithmetic, only : ieee_is_nan
  use mo_mpr_util_types, only : mpr_get_coord_alias, MprBaseType, mpr_add_coord_alias
  use mo_mpr_utils, only : get_n_initialized, get_index_in_vector, read_coordinate_from_file, read_variable_from_file
  use mo_mpr_file, only : filenameNamelistMprDefault
  use mo_append, only : append
  use flogging

  implicit none

  private

  public :: get_index_in_coordinate, add_coordinate, create_upscaler_name, create_upscaler_alias, newCoordinate, &
            combine_coordinates, create_integral_coordinate, split_coordinate, get_common_unit, check_within_bounds

  ! --------------------------------------------------------------------------------------

  ! TODO: might be moved to other place and made available to data arrays also
  !> type containing all information of attributes that can be attached to any (coordinate) variable
  type :: Attribute
    private
    ! name of attribute
    character(maxNameLength) :: name
    ! type, see https://github.com/Unidata/netcdf-fortran/blob/9597d59f14ab13f377e43733fc65f74e51aca705/fortran/netcdf_constants.f90
    integer(i4) :: type
    ! length of characters or of vector if 1d
    integer(i4) :: length
    ! value containers for all type and scalar/1d combinations possible
    character(maxStringLength) :: valuechar
    integer(i1) :: valuei1
    integer(i2) :: valuei2
    integer(i4) :: valuei4
    real(sp) :: valuesp
    real(dp) :: valuedp
    integer(i1), dimension(:), allocatable :: valuesi1
    integer(i2), dimension(:), allocatable :: valuesi2
    integer(i4), dimension(:), allocatable :: valuesi4
    real(sp), dimension(:), allocatable :: valuessp
    real(dp), dimension(:), allocatable :: valuesdp
  contains
    procedure :: get_valuechar
    procedure :: initchar
    procedure :: initi1
    procedure :: initi2
    procedure :: initi4
    procedure :: initsp
    procedure :: initdp
    procedure :: initsi1
    procedure :: initsi2
    procedure :: initsi4
    procedure :: initssp
    procedure :: initsdp
    generic, public :: get_value => get_valuechar
    generic, public :: init => initchar, initi1, initi2, initi4, initsp, initdp, &
            initsi1, initsi2, initsi4, initssp, initsdp
    procedure, public :: get_type
    procedure, public :: get_length
    procedure, public :: get_name
  end type Attribute


  type, extends(MprBaseType), public :: Coordinate

    real(dp) :: step
    integer(i8) :: count
    integer(i4) :: rank
    integer(i4) :: corners
    real(dp), dimension(:), allocatable :: values
    real(dp), dimension(:), allocatable :: centersCoord1
    real(dp), dimension(:), allocatable :: centersCoord2
    real(dp), dimension(:,:), allocatable :: cornersCoord1
    real(dp), dimension(:,:), allocatable :: cornersCoord2
    integer(i4), dimension(:), allocatable :: nodes
    integer(i4) :: staggerId
    character(6) :: staggerName
    type(Attribute), dimension(maxNoAttributes) :: attributes
    real(dp), dimension(2) :: bounds
    character(maxStringLength):: projString
    character(maxStringLength):: unit
    character(maxNameLength), dimension(:), allocatable :: subDims
    integer(i8), dimension(:), allocatable :: subDimSizes

  contains
    private
    procedure, public :: is_finalized
    procedure, public :: reset
    procedure :: fromFile
    procedure, public :: from_other
    procedure :: from_range
    procedure :: from_values
    procedure :: set_values
    procedure :: check_bound
    procedure :: set_step
    procedure :: set_count
    procedure :: is_equal_coordinate
    procedure :: is_unequal_coordinate
    procedure :: insert_bound
    procedure :: get_bound_index
    procedure :: set_bound
    procedure :: get_attributes_from_nc_var
    procedure, public :: is_ascending
    procedure, public :: set_2d_count
    procedure, public :: set_polygons_from_2d
    procedure, public :: is_2d
    procedure, public :: is_polygon
    procedure, public :: write_coordinate
    procedure, public :: get_bound
    procedure, public :: get_stats
    procedure, public :: has_attribute
    procedure, public :: get_attribute
    generic, public :: operator(.eq.) => is_equal_coordinate
    generic, public :: operator(.ne.) => is_unequal_coordinate
  end type Coordinate

  interface Coordinate
    procedure newCoordinate
  end interface Coordinate

  type, public :: CoordinatePointer
    ! this is a pointer to a Coordinate
    type(Coordinate), pointer :: coord_p => null()
  contains
    private
    procedure :: set_coordinate_pointer_by_name
    procedure :: set_coordinate_pointer_by_id
    procedure :: set_coordinate_pointer_by_dim
    generic, public :: set_coordinate_pointer => &
            set_coordinate_pointer_by_name, &
            set_coordinate_pointer_by_id, &
            set_coordinate_pointer_by_dim

  end type CoordinatePointer

  ! each Coordinate type is allocatable
  type(Coordinate), dimension(MaxNoCoords), public, target :: MPR_COORDINATES ! all coordinates

  interface create_upscaler_name
    procedure create_upscaler_name_dim, create_upscaler_name_dimp
  end interface create_upscaler_name

  interface create_upscaler_alias
    procedure create_upscaler_alias_dim, create_upscaler_alias_dimp
  end interface create_upscaler_alias

  ! --------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------
contains

  subroutine initchar(att, name, arg)
    !< init Attribute with character value
    character(*), intent(in) :: name
    character(*), intent(in) :: arg
    class(Attribute), intent(out) :: att
    
    att%name = name
    att%type = 2_i4
    att%length = len_trim(arg)
    att%valuechar = arg
    
  end subroutine initchar
  
  subroutine initi1(att, name, arg)
    !< init Attribute with byte value
    character(*), intent(in) :: name
    integer(i1), intent(in) :: arg
    class(Attribute), intent(out) :: att
    
    att%name = name
    att%type = 1_i4
    att%length = 1_i4
    att%valuei1 = arg
    
  end subroutine initi1
  
  subroutine initi2(att, name, arg)
    !< init Attribute with integer (i2) value
    character(*), intent(in) :: name
    integer(i2), intent(in) :: arg
    class(Attribute), intent(out) :: att
    
    att%name = name
    att%type = 3_i4
    att%length = 1_i4
    att%valuei2 = arg
    
  end subroutine initi2
  
  subroutine initi4(att, name, arg)
    !< init Attribute with integer (i4) value
    character(*), intent(in) :: name
    integer(i4), intent(in) :: arg
    class(Attribute), intent(out) :: att
    
    att%name = name
    att%type = 4_i4
    att%length = 1_i4
    att%valuei4 = arg
    
  end subroutine initi4
  
  subroutine initsp(att, name, arg)
    !< init Attribute with short (sp) value
    character(*), intent(in) :: name
    real(sp), intent(in) :: arg
    class(Attribute), intent(out) :: att
    
    att%name = name
    att%type = 5_i4
    att%length = 1_i4
    att%valuesp = arg
    
  end subroutine initsp
  
  subroutine initdp(att, name, arg)
    !< init Attribute with double (sp) value
    character(*), intent(in) :: name
    real(dp), intent(in) :: arg
    class(Attribute), intent(out) :: att
    
    att%name = name
    att%type = 6_i4
    att%length = 1_i4
    att%valuedp = arg
    
  end subroutine initdp
  
  subroutine initsi1(att, name, arg)
    !< init Attribute with short byte vector
    character(*), intent(in) :: name
    integer(i1), intent(in), dimension(:) :: arg
    class(Attribute), intent(out) :: att
    
    att%name = name
    att%type = 1_i4
    att%length = size(arg)
    att%valuesi1 = arg
    
  end subroutine initsi1
  
  subroutine initsi2(att, name, arg)
    !< init Attribute with short i2 vector
    character(*), intent(in) :: name
    integer(i2), intent(in), dimension(:) :: arg
    class(Attribute), intent(out) :: att
    
    att%name = name
    att%type = 3_i4
    att%length = size(arg)
    att%valuesi2 = arg
    
  end subroutine initsi2
  
  subroutine initsi4(att, name, arg)
    !< init Attribute with short i4 vector
    character(*), intent(in) :: name
    integer(i4), intent(in), dimension(:) :: arg
    class(Attribute), intent(out) :: att
    
    att%name = name
    att%type = 4_i4
    att%length = size(arg)
    att%valuesi4 = arg
    
  end subroutine initsi4
  
  subroutine initssp(att, name, arg)
    !< init Attribute with short sp vector
    character(*), intent(in) :: name
    real(sp), intent(in), dimension(:) :: arg
    class(Attribute), intent(out) :: att
    
    att%name = name
    att%type = 5_i4
    att%length = size(arg)
    att%valuessp = arg
    
  end subroutine initssp
  
  subroutine initsdp(att, name, arg)
    !< init Attribute with short dp vector
    character(*), intent(in) :: name
    real(dp), intent(in), dimension(:) :: arg
    class(Attribute), intent(out) :: att

    att%name = name
    att%type = 6_i4
    att%length = size(arg)
    att%valuesdp = arg
    
  end subroutine initsdp
  
  function get_type(att) result(type)
    !< get type of attribute
    !< 1-byte, 2-char, 3-short, 4-int, 5-float, 6-double
    class(Attribute), intent(in) :: att
    integer(i4) :: type

    type = att%type

  end function get_type

  function get_length(att) result(length)
    !< get length of attribute
    class(Attribute), intent(in) :: att
    integer(i4) :: length

    length = att%length

  end function get_length

  function get_name(att) result(name)
    !< get name of attribute
    class(Attribute), intent(in) :: att
    character(maxNameLength) :: name

    name = att%name

  end function get_name

  function get_valuechar(att) result(value)
    !< get character value of attribute
    class(Attribute), intent(in) :: att
    character(maxNameLength) :: value

    value = ''
    if (att%type == 2_i4) value = att%valuechar

  end function get_valuechar

  subroutine set_coordinate_pointer_by_name(self, dimIdentifier, ierr)
    class(CoordinatePointer), intent(inout) :: self
    character(*), intent(in) :: dimIdentifier
    integer(i4), intent(out), optional :: ierr
    integer(i4), dimension(1) :: coordIndex
    integer(i4) :: nCoords

    if (present(ierr)) ierr = 0_i4
    coordIndex = get_index_in_coordinate(dimIdentifier)
    nCoords = get_n_initialized(MPR_COORDINATES)
    if (coordIndex(1) > nCoords) then
      if (present(ierr)) then
        ierr = 1_i4
      else
        log_error(*) "Coordinate name '", trim(dimIdentifier), "' is not defined explicitly.", &
                " Add the target coordinate or an appropriate alias in the file ", trim(filenameNamelistMprDefault)
        stop 1
      end if
      return
    end if

    self%coord_p => MPR_COORDINATES(coordIndex(1))

  end subroutine set_coordinate_pointer_by_name

  subroutine set_coordinate_pointer_by_id(self, dimIdentifier, ierr)
    class(CoordinatePointer), intent(inout) :: self
    integer(i4), intent(in) :: dimIdentifier
    integer(i4), intent(out), optional :: ierr
    integer(i4) :: nCoords
    character(maxStringLength) :: errorMessage

    if (present(ierr)) ierr = 0_i4
    nCoords = get_n_initialized(MPR_COORDINATES)
    if (dimIdentifier > nCoords) then
      errorMessage = "Coordinate id '"//trim(compress(num2str(dimIdentifier)))//"' is not found in target coordinates."// &
              " Add the target coordinate or an appropriate alias in the file "// trim(filenameNamelistMprDefault)
      if (present(ierr)) then
        ierr = 1_i4
      else
        log_error(*) trim(errorMessage)
        stop 1
      end if
      return
    end if

    self%coord_p => MPR_COORDINATES(dimIdentifier)

  end subroutine set_coordinate_pointer_by_id

  subroutine set_coordinate_pointer_by_dim(self, dimIdentifier, ierr)
    class(CoordinatePointer), intent(inout) :: self
    type(Coordinate), intent(in) :: dimIdentifier
    integer(i4), intent(out), optional :: ierr
    integer(i4) :: nCoords, iCoord, foundCoord
    character(maxStringLength) :: errorMessage

    if (present(ierr)) ierr = 0_i4
    foundCoord = 0_i4
    nCoords = get_n_initialized(MPR_COORDINATES)
    do iCoord = 1, nCoords
      if (dimIdentifier == MPR_COORDINATES(iCoord)) then
        foundCoord = iCoord
      end if
    end do
    if (foundCoord == 0_i4) then
      errorMessage = "Coordinate dim '"//trim(dimIdentifier%name)//"' is not found in target coordinates."// &
              " Add the target coordinate or an appropriate alias in the file "//trim(filenameNamelistMprDefault)
      if (present(ierr)) then
        ierr = 1_i4
      else
        log_error(*) trim(errorMessage)
        stop 1
      end if
      return
    end if

    self%coord_p => MPR_COORDINATES(foundCoord)

  end subroutine set_coordinate_pointer_by_dim

  function newCoordinate(name, id, stagger, fileName, values, bound, start, step, count_, &
          attributeNames, attributeValuesChar, projString, subDims, unit, &
          centersCoord1, centersCoord2, cornersCoord1, cornersCoord2)
    character(*), intent(in) :: name
    integer(i4), intent(in) :: id
    character(*), intent(in) :: stagger
    character(*), intent(in), optional :: fileName
    real(dp), dimension(:), intent(in), optional :: values
    real(dp), intent(in), optional :: bound
    real(dp), intent(in), optional :: start
    real(dp), intent(in), optional :: step
    integer(i8), intent(in), optional :: count_
    character(*), dimension(:), intent(in), optional :: attributeNames
    character(*), dimension(:), intent(in), optional :: attributeValuesChar
    character(*), intent(in), optional :: projString
    character(*), dimension(:), intent(in), optional :: subDims
    character(*), intent(in), optional :: unit
    real(dp), dimension(:), intent(in), optional :: centersCoord1
    real(dp), dimension(:), intent(in), optional :: centersCoord2
    real(dp), dimension(:,:), intent(in), optional :: cornersCoord1
    real(dp), dimension(:,:), intent(in), optional :: cornersCoord2
    type(Coordinate) :: newCoordinate

    character(maxStringLength) :: fileNameArg
    real(dp), dimension(:), allocatable :: valuesArg
    real(dp) :: boundArg
    real(dp) :: startArg
    real(dp) :: stepArg
    integer(i4) :: iPolygon, iAtt
    integer(i8) :: countArg
    logical, dimension(:,:), allocatable :: mask

    ! check values and set
    if (trim(name) == trim(defaultAlias)) then
      log_error(*) "Error during initialization. The coordinate has no valid name."
      stop 1
    else
      newCoordinate%name = trim(name)
      newCoordinate%id = id
      newCoordinate%is_initialized = .true.
    end if

    newCoordinate%staggerName = trim(stagger)
    select case(trim(stagger))
    case("start")
      newCoordinate%staggerId = 0_i8
    case("center")
      newCoordinate%staggerId = 1_i8
    case("end")
      newCoordinate%staggerId = 2_i8
    case default
      log_error(*) "Error during initialization of Coordinate '", trim(name), "'. The stagger string '", &
              trim(stagger),"' is not known. Use 'start', 'center' or 'end'."
      stop 1
    end select

    ! init values (all those can be nodata_dp...)
    ! init all defaults
    fileNameArg = defaultAlias
    if (present(fileName)) then
      fileNameArg = fileName
    end if
    if (present(values)) then
      valuesArg = values
    else
      allocate(valuesArg(0))
    end if
    stepArg = nodata_dp
    if (present(step)) then
      stepArg = step
    end if
    startArg = nodata_dp
    if (present(start)) then
      startArg = start
    end if
    countArg = defaultCoordCount
    if (present(count_)) then
      countArg = count_
    end if
    boundArg = nodata_dp
    if (present(bound)) then
      boundArg = bound
    end if
    newCoordinate%step = nodata_dp
    newCoordinate%bounds = nodata_dp
    newCoordinate%count = defaultCoordCount
    newCoordinate%projString = defaultAlias
    newCoordinate%unit = defaultCoordUnits
    newCoordinate%rank = 1_i4
    newCoordinate%corners = 2_i4

    if (present(unit)) then
      newCoordinate%unit = unit
    end if

    ! set the properties
    ! first path: init from file
    if (trim(fileNameArg) /= trim(defaultAlias)) then
      call newCoordinate%fromFile(fileNameArg, boundArg)
    ! second path: init from values (but only if valid - at least one non-Nan value))
    else if (size(valuesArg) > 0_i4) then
      if (.not. all(eq(valuesArg, nodata_dp))) then
        call newCoordinate%from_values(valuesArg, boundArg)
      ! third path: init from range (but only if valid)
      else if (ne(stepArg, nodata_dp)) then
        if (ne(startArg, nodata_dp) .and. countArg /= defaultCoordCount) then
          call newCoordinate%from_range(stepArg, startArg, countArg)
        ! fourth path: set only single attributes for init from other later on
        else
          call newCoordinate%set_step(step=stepArg)
        end if
      ! fourth path: set only single attributes for init from other later on
      else if (ne(boundArg, nodata_dp)) then
        call newCoordinate%set_bound(boundArg)
      ! fifth path: set only single attributes for init from other later on
      else if (countArg /= defaultCoordCount) then
        call newCoordinate%set_count(countArg)
      end if
    ! third path: init from range (but only if valid)
    else if (ne(stepArg, nodata_dp)) then
      if (ne(startArg, nodata_dp) .and. countArg /= defaultCoordCount) then
        call newCoordinate%from_range(stepArg, startArg, countArg)
      ! fourth path: set only single attributes for init from other later on
      else
        call newCoordinate%set_step(step=stepArg)
      end if
    ! fourth path: set only single attributes for init from other later on
    else if (ne(boundArg, nodata_dp)) then
      call newCoordinate%set_bound(boundArg)
    ! fifth path: set only single attributes for init from other later on
    else if (countArg /= defaultCoordCount) then
      call newCoordinate%set_count(countArg)
    end if

    ! set the attributes
    if (present(attributeNames) .and. present(attributeValuesChar)) then
      if (size(attributeNames) == size(attributeValuesChar)) then
        if (size(attributeNames) > maxNoAttributes) then
          log_error(("(1X,A,I0,A,I0,A)")) 'Provided more attributes (', size(attributeNames) , &
              ')than allowed (', maxNoAttributes, ').'
          stop 1
        end if
        do iAtt=1, size(attributeNames)
          call newCoordinate%attributes(iAtt)%init(attributeNames(iAtt), attributeValuesChar(iAtt))
        end do
      end if
    end if

    ! set the projection information
    if (present(projString)) then
      newCoordinate%projString = projString
    end if

    ! set the subDims information
    if (present(subDims)) then
      log_debug(*) 'newCoordinate: for Coordinate ', trim(newCoordinate%name), ' initializing subDims: ', subDims
      if (size(subDims) > 0) then
        newCoordinate%subDims = subDims
        newCoordinate%rank = size(subDims)
        if (newCoordinate%is_finalized() .and. newCoordinate%count == defaultCoordCount) then
          call newCoordinate%set_2d_count()
        end if
      end if
    end if
    if (.not. allocated(newCoordinate%subDimSizes)) allocate(newCoordinate%subDimSizes(0_i8))
    if (.not. allocated(newCoordinate%subDims)) allocate(newCoordinate%subDims(0_i4))

    ! store the polygon information
    if (present(centersCoord1) .and. present(centersCoord2)) then
      newCoordinate%centersCoord1 = centersCoord1
      newCoordinate%centersCoord2 = centersCoord2
      call newCoordinate%set_count(size(centersCoord2, kind=i8))
      newCoordinate%rank = 2_i4
      if (present(cornersCoord1) .and. present(cornersCoord2)) then
        newCoordinate%cornersCoord1 = cornersCoord1
        newCoordinate%cornersCoord2 = cornersCoord2
        newCoordinate%corners = size(cornersCoord2,dim=1)
        allocate(newCoordinate%nodes(size(cornersCoord2,dim=1)))
        mask = ne(cornersCoord2, nodata_dp)
        if (.not. all(mask)) then
          do iPolygon=1, size(cornersCoord2, dim=2)
            newCoordinate%nodes(iPolygon) = count(mask(:, iPolygon))
          end do
        else
          newCoordinate%nodes = size(cornersCoord2,dim=1)
        end if
      end if
    end if

    if (.not. allocated(newCoordinate%nodes)) allocate(newCoordinate%nodes(0_i4))

    ! if not properly initialized, issue warning
    if (.not. newCoordinate%is_finalized()) then
      log_debug(*) "Cannot properly initialize coordinate '", trim(newCoordinate%name), "'."
    ! else if (.not. newCoordinate%is_polygon()) then
    !   allocate(newCoordinate%centersCoord1(newCoordinate%count))
    !   allocate(newCoordinate%centersCoord2(newCoordinate%count))
    !   allocate(newCoordinate%cornersCoord1(newCoordinate%corners, newCoordinate%count))
    !   allocate(newCoordinate%cornersCoord2(newCoordinate%corners, newCoordinate%count))
    end if


  end function newCoordinate

  logical function is_equal_coordinate(coord1, coord2)
    !< checks wheter two Coordinates are equal
    !< this is the case if
    !< 1) the names are equal or in the same CoordAlias group AND
    !< 2) the stagger is equal AND
    !< 3) both instances are finalized and the values of the vector are equal within a tolerance

    !> reference coordinate
    class(Coordinate), intent(in) :: coord1
    !> coordinate to compare
    class(Coordinate), intent(in) :: coord2

    character(maxNameLength), dimension(:), allocatable :: aliases
    logical, dimension(:), allocatable :: isNotAliases
    integer(i4) :: i, iSubCoord

    ! let's be negative
    is_equal_coordinate = .false.

    ! name should be equal...
    if (coord1%name /= coord2%name) then
      ! ... or same alias
      ! look in aliases LUTs if found there
      call mpr_get_coord_alias(coord2%name, .false., aliases)
      if (size(aliases) == 0_i4) then
        ! coord2%name is not even in any LUT
        log_subtrace(*) 'equalCoordinate: different names and no alias'
        return
      else
        allocate(isNotAliases(size(aliases)))
        isNotAliases = (/(aliases(i) /= coord1%name, i = 1, size(aliases))/)
        if (all(isNotAliases)) then
          ! coord1%name is not in coord2%name's LUT
          log_subtrace(*) 'equalCoordinate: different names and alias groups do not match'
          return
        end if
      end if
    end if

    ! data reference should also be equal
    if (coord1%staggerId /= coord2%staggerId) then
      log_subtrace(*) 'equalCoordinate: different stagger'
      return
    end if

    ! compare the actual values (needs to be finalized)
    if (coord1%is_finalized() .and. coord2%is_finalized()) then
      ! check for type
      if (coord1%is_2d() .and. coord2%is_2d()) then
        if (coord1%is_polygon() .and. coord2%is_polygon()) then
          ! ...comparing two polygon 2d coordinates,
          if (coord1%corners /= coord2%corners) then
            ! different size of values
            log_subtrace(*) 'equalCoordinate: different corners size'
            return
          end if
          if (allocated(coord1%cornersCoord1) .and. allocated(coord2%cornersCoord1)) then
            if (any(abs(coord1%cornersCoord1 - coord2%cornersCoord1) > maxTolerance) .or. &
                    any(abs(coord1%cornersCoord2 - coord2%cornersCoord2) > maxTolerance)) then
              ! different values on values
              log_subtrace(*) 'equalCoordinate: different corners values'
              return
            end if
          end if
        else if (.not. coord1%is_polygon() .and. .not. coord2%is_polygon()) then
          ! ...comparing two regular 2d coordinates,
          if (size(coord1%subDims) /= size(coord2%subDims)) then
            ! different size of subDims
            log_subtrace(*) 'equalCoordinate: 2d coordinate different number of subDims'
            return
          end if
          if (any([(coord1%subDims(iSubCoord) /= coord2%subDims(iSubCoord), iSubCoord=1, size(coord1%subDims))])) then
            ! different subDims
            log_subtrace(*) 'equalCoordinate: 2d coordinate subDims are not equal'
            return
          end if

        else
          ! two different types of 2d coordinates
          return
        end if
      else if (.not. coord1%is_2d() .and. .not. coord2%is_2d()) then
        ! ...comparing two 1d coordinates, data values should also be equal
        if (size(coord1%values, kind=i8) /= size(coord2%values, kind=i8)) then
          ! different size of values
          log_subtrace(*) 'equalCoordinate: different vector size'
          return
        end if
        if (any(abs(coord1%values - coord2%values) > maxTolerance)) then
          ! different values on values
          log_subtrace(*) 'equalCoordinate: different values, max deviation: ', maxval(abs(coord1%values - coord2%values))
          return
        end if
      else
        ! comparing 1d and 2d coordinate coordinates
        log_subtrace(*) 'equalCoordinate: comparing finalized 1d and 2d coordinates'
        return
      end if
    else
      log_subtrace(*) 'equalCoordinate: at least one coordinate is not finalized'
      return
    end if
    is_equal_coordinate = .true.

  end function is_equal_coordinate

  logical function is_unequal_coordinate(coord1, coord2)
    !< checks wheter two Coordinates are unequal

    !> reference coordinate
    class(Coordinate), intent(in) :: coord1
    !> coordinate to compare
    class(Coordinate), intent(in) :: coord2

    is_unequal_coordinate = (.not. coord1 == coord2)

  end function is_unequal_coordinate

  logical function is_2d(dim)
    class(Coordinate), intent(in) :: dim

    is_2d = .false.
    if (dim%rank > 1) then
      is_2d = .true.
    end if

  end function is_2d

  logical function is_polygon(dim)
    class(Coordinate), intent(in) :: dim

    is_polygon = .false.

    if (dim%is_2d() .and. dim%corners > 2) then
      is_polygon = .true.
    end if

  end function is_polygon

  subroutine fromFile(self, fileName, boundArg)
    !< initialize a Coordinate instance based on an existing Coordinate in a file
    !< currently this works for a coordinate variable if it is 1D or 2D on regular, non-shifted grid

    !> Coordinate type-bound procedure
    class(Coordinate), intent(inout) :: self
    !> the file name of the netcdf file containing the coordinate self%name
    character(*), intent(in) :: fileName
    !> bound to be set (optional - needed only, when there is an irregular step size)
    real(dp), intent(in), optional :: boundArg

    type(NcDataset) :: nc
    type(NcVariable) :: ncVar
    type(NcDimension) :: nc_dim
    type(NcDimension), dimension(:), allocatable :: NcDims
    integer(i4), dimension(:), allocatable :: varShape
    real(dp), dimension(:, :), allocatable :: data2D, dummy
    character(maxNameLength), dimension(:), allocatable :: aliases, coordNames
    character(maxNameLength) :: ncTitle
    real(dp), dimension(:), allocatable :: tempValues
    real(dp) :: bound
    logical, dimension(:, :), allocatable :: ncMask

    integer(i4) :: j, vectorSize, coordIndex, dimSize
    integer(i4), allocatable, dimension(:) :: dimSizes
    integer(i8) :: i
    character(maxNameLength) :: errorString, boundVariable

    ! init
    bound = nodata_dp
    if (present(boundArg)) then
      if (.not. ieee_is_nan(boundArg)) bound = boundArg
    end if

    log_debug(*) "Initializing coordinate '", trim(self%name), "' from file."

    ! read Dataset and store the variable of interest
    nc = NcDataset(fileName, "r")
    if (nc%hasAttribute('title')) then
      ! SCRIP format
      call nc%getAttribute('title', ncTitle)
      self%name = trim(ncTitle)
      ! found the title, now read all the different field information
      ! start with the coordinates
      call read_coordinate_from_file(nc, 'grid_size', dimSize)
      call self%set_count(int(dimSize, kind=i8))
      ! assuming a rank of 2
      ! SCRIP uses rank 1 for unstructured grids, we use it differently
      ! call read_coordinate_from_file(nc, 'grid_rank', self%rank)
      self%rank = 2_i4
      self%subDims =['lon', 'lat']
      call read_coordinate_from_file(nc, 'grid_corners', self%corners)

      allocate(self%cornersCoord1(self%corners, self%count))
      allocate(self%cornersCoord2(self%corners, self%count))
      allocate(self%centersCoord1(self%count))
      allocate(self%centersCoord2(self%count))
      allocate(self%nodes(self%count))

      call read_variable_from_file(nc, 'grid_center_lon', self%centersCoord1)
      call read_variable_from_file(nc, 'grid_center_lat', self%centersCoord2)
      call read_variable_from_file(nc, 'grid_corner_lon', self%cornersCoord1)
      ! reading also the mask to infer the number of corners/nodes per grid cell
      call read_variable_from_file(nc, 'grid_corner_lat', self%cornersCoord2, ncMask)
      if (.not. all(ncMask)) then
        do i = 1, self%count
          self%nodes(i) = self%corners - count(.not. ncMask(:, i))
        end do
      else
        ! set the number of nodes to the maximum number of corners
        self%nodes = self%corners
      end if
      deallocate(ncMask)

      ! for now ignore the mask as this applied later on with the data array only
      ! call self%read_mask_from_file(nc, 'grid_imask', self%mask)
      call read_variable_from_file(nc, 'grid_dims', dimSizes)
      self%subDimSizes = int(dimSizes, kind=i8)

      ncVar = nc%getVariable('grid_center_lon')
      ! this is mimicking the retrieval of grids.f in SCRIP -> assume unit to be similar over all variables
      call ncVar%getAttribute('units', self%unit)

      log_subtrace(*) 'Reading from SCRIP dimension:'
      log_subtrace(*) 'grid_size ', self%count
      log_subtrace(*) 'grid_rank ', self%rank
      log_subtrace(*) 'grid_corners ', self%corners
      log_subtrace(*) 'grid_corner_lon ', self%cornersCoord1
      log_subtrace(*) 'grid_corner_lat ', self%cornersCoord2
      log_subtrace(*) 'grid_center_lon ', self%centersCoord1
      log_subtrace(*) 'grid_center_lat ', self%centersCoord2
      log_subtrace(*) 'grid_dims ', self%subDimSizes
      log_subtrace(*) 'unit ', trim(self%unit)
    else if (nc%hasVariable(trim(self%name))) then
      ! the coordinate is actually a variable in the NcDataset
      ncVar = nc%getVariable(trim(self%name))
      varShape = ncVar%getShape()
      ! infer the dimensions...
      select case(size(varShape))
      case(1)
        ! the variable is dependent on 1 dimension only
        allocate(tempValues(varShape(1)))
        call ncVar%getData(tempValues)

        if (ncVar%hasAttribute('bounds')) then
          call ncVar%getAttribute('bounds', boundVariable)
          ncVar = nc%getVariable(trim(boundVariable))
          call ncVar%getData(dummy)
          if (self%staggerId >= 1_i8) then
            ! end- or center-bound values
            bound = dummy(1, 1)
            ! overwrite data by bound information
            ! tempValues = dummy(2, 1:size(dummy,2))
          else
            bound = dummy(2, size(dummy,2))
            ! overwrite data by bound information
            ! tempValues = dummy(1, 1:size(dummy,2))
          end if
          deallocate(dummy)
          log_trace(*) "fromFile: read values from auxiliary bound coordinate"
          log_trace(*) "fromFile: bound: ", bound
        else
          log_trace(*) "fromFile: read values from coordinate"
        end if
        log_subtrace(*) "fromFile: stagger: ", trim(self%staggerName)
        log_subtrace(*) "fromFile: tempValues: ", tempValues
      case(2)
        ! the variable is dependent on 2 dimensions
        NcDims = ncVar%getDimensions()
        allocate(coordNames(size(NcDims)))
        do i = 1, size(NcDims)
          coordNames(i) = trim(NcDims(i)%getName())
        end do

        call mpr_get_coord_alias(self%name, .true., aliases)
        if (size(aliases) == 0_i4) then
          errorString = ""
          do i = 1, size(coordNames)
            errorString = trim(errorString) // trim(coordNames(i))
          end do
          log_error(*)  "The coordinate variable ", trim(self%name), " has no coordinate alias. Add one of '", &
                  errorString, "' to the variables' aliases in the file ", &
                  filenameNamelistMprDefault
          call nc%close()
          stop 1
        end if

        ! check for all aliases if it is a coordinate, get it
        coordIndex = 0_i4
        alias : do i = 1, size(aliases)
          do j = 1, size(coordNames)
            if (coordNames(j) == aliases(i)) then
              nc_dim = nc%getDimension(aliases(i))
              vectorSize = nc_dim%getLength()
              coordIndex = j
              exit alias
            end if
          end do
        end do alias
        deallocate(coordNames)
        if (coordIndex == 0_i4) then
          log_error(*)  "The coordinate variable ", trim(self%name), " or any of its aliases is no coordinate ", &
                  "in the nc file ", trim(fileName)
          call nc%close()
          stop 1
        end if
        allocate(tempValues(vectorSize))
        ! now check if the coordinate is regular (grid) and set the vector values
        call ncVar%getData(data2D)
        if (coordIndex == 1_i4) then
          do i = 1, size(tempValues)
            if (all(eq(data2D(i, :), data2D(i, 1)))) then
              tempValues(i) = data2D(i, 1)
            else
              log_error(*)  "The coordinate variable ", trim(self%name), " is not orthogonal but tilted ", &
                      "in the nc file ", trim(fileName)
              call nc%close()
              stop 1
            end if
          end do
        else
          do i = 1, size(tempValues)
            if (all(eq(data2D(:, i), data2D(1, i)))) then
              tempValues(i) = data2D(1, i)
            else
              log_error(*)  "The coordinate variable ", trim(self%name), " is not orthogonal but tilted ", &
                      "in the nc file ", trim(fileName)
              call nc%close()
              stop 1
            end if
          end do
        end if
      case default
        log_error("(1X,A,I0,A)") size(varShape), "-dimensional coordinate arrays are not supported as coordinates."
        call nc%close()
        stop 1
      end select
      ! get all the attribute names of the original coordinate variable (change, might be set to bnds)
      ncVar = nc%getVariable(trim(self%name))
      call self%get_attributes_from_nc_var(ncVar)
      ! set special attribute "units"
      if (ncVar%hasAttribute('units')) then
        call ncVar%getAttribute('units', self%unit)
      end if
    else
      log_debug(*) "Cannot properly initialize coordinate '", trim(self%name), "'."
      call nc%close()
      return
    end if

    call nc%close()

    if (.not. self%is_2d()) then
      ! set step not considering the bound
      call self%set_step(values=tempValues)
      if (ne(bound, nodata_dp)) then
        ! everything is provided, init the values
        call self%set_values(tempValues, bound)
        ! set step considering the bound
        call self%set_step(values=self%values(0:self%count))
      else
        ! try to init the values without a bound
        call self%set_values(tempValues)
      end if
    end if

  end subroutine fromFile

  subroutine get_attributes_from_nc_var(self, ncVar)
    !< get all the attributes from a netcdf variable and store them in the self%attributes property
    class(Coordinate), intent(inout) :: self
    type(ncVariable), intent(in) :: ncVar

    character(maxNameLength), dimension(:), allocatable :: attributeNames
    integer(i4) :: iAtt, attType, attLength

    character(maxStringLength) :: valuechar
    integer(i1) :: valuei1
    integer(i2) :: valuei2
    integer(i4) :: valuei4
    real(sp) :: valuesp
    real(dp) :: valuedp
    integer(i1), dimension(:), allocatable :: valuesi1
    integer(i2), dimension(:), allocatable :: valuesi2
    integer(i4), dimension(:), allocatable :: valuesi4
    real(sp), dimension(:), allocatable :: valuessp
    real(dp), dimension(:), allocatable :: valuesdp
    logical :: hasAttr

    ! get all the names
    attributeNames = ncVar%getAttributeNames()

    do iAtt=1, size(attributeNames)
      ! get the type and length
      hasAttr = ncVar%hasAttribute(attributeNames(iAtt), attType, attLength)
      ! check for 1d vector values being longer than allowed
      if (hasAttr .and. attLength > maxNoAttributeValues .and. attType /= 2_i4) then
        log_error("(1X,A,I0,A,A,A,A)") 'Received more than ', maxNoAttributeValues, &
                ' attribute values for attribute ', trim(attributeNames(iAtt)), &
                ' of variable', trim(ncVar%getName())
        stop 1
      end if

      select case(attType)
      ! 1-byte, 2-char, 3-short, 4-int, 5-float, 6-double
      case(1_i4)
        if (attLength == 1) then
          call ncVar%getAttribute(attributeNames(iAtt), valuei1)
          call self%attributes(iAtt)%init(attributeNames(iAtt), valuei1)
        else
          call ncVar%getAttribute(attributeNames(iAtt), valuesi1)
          call self%attributes(iAtt)%init(attributeNames(iAtt), valuesi1)
        end if
      case(2_i4)
        call ncVar%getAttribute(attributeNames(iAtt), valuechar)
        call self%attributes(iAtt)%init(attributeNames(iAtt), valuechar)
      case(3_i4)
        if (attLength == 1) then
          call ncVar%getAttribute(attributeNames(iAtt), valuei2)
          call self%attributes(iAtt)%init(attributeNames(iAtt), valuei2)
        else
          call ncVar%getAttribute(attributeNames(iAtt), valuesi2)
          call self%attributes(iAtt)%init(attributeNames(iAtt), valuesi2)
        end if
      case(4_i4)
        if (attLength == 1) then
          call ncVar%getAttribute(attributeNames(iAtt), valuei4)
          call self%attributes(iAtt)%init(attributeNames(iAtt), valuei4)
        else
          call ncVar%getAttribute(attributeNames(iAtt), valuesi4)
          call self%attributes(iAtt)%init(attributeNames(iAtt), valuesi4)
        end if
      case(5_i4)
        if (attLength == 1) then
          call ncVar%getAttribute(attributeNames(iAtt), valuesp)
          call self%attributes(iAtt)%init(attributeNames(iAtt), valuesp)
        else
          call ncVar%getAttribute(attributeNames(iAtt), valuessp)
          call self%attributes(iAtt)%init(attributeNames(iAtt), valuessp)
        end if
      case(6_i4)
        if (attLength == 1) then
          call ncVar%getAttribute(attributeNames(iAtt), valuedp)
          call self%attributes(iAtt)%init(attributeNames(iAtt), valuedp)
        else
          call ncVar%getAttribute(attributeNames(iAtt), valuesdp)
          call self%attributes(iAtt)%init(attributeNames(iAtt), valuesdp)
        end if
      end select
    end do

  end subroutine get_attributes_from_nc_var
  
  subroutine from_values(self, values, boundArg)
    !< initialize a Coordinate instance based on a vector of values

    !> Coordinate type-bound procedure
    class(Coordinate), intent(inout) :: self
    !> vector of values to be set
    real(dp), dimension(:), intent(in) :: values
    !> bound to be set (optional - needed only, when there is an irregular step size)
    real(dp), intent(in), optional :: boundArg

    real(dp) :: bound

    ! init
    bound = nodata_dp
    if (present(boundArg)) then
      bound = boundArg
    end if

    log_debug(*) "Initializing coordinate '", trim(self%name), "' from values."

    ! set step not considering the bound
    call self%set_step(values=values)
    if (ne(bound, nodata_dp)) then
      ! everything is provided, init the values
      call self%set_values(values, bound)
      ! set step considering the bound
      call self%set_step(values=self%values)
    else
      ! try to init the values without a bound
      call self%set_values(values)
    end if

  end subroutine from_values

  subroutine from_range(self, step, start, count)
    !< initialize a Coordinate instance based on a defined range

    !> Coordinate type-bound procedure
    class(Coordinate), intent(inout) :: self
    !> step or increment of value range
    real(dp), intent(in) :: step
    !> start value
    real(dp), intent(in):: start
    !> number of values in range
    integer(i8), intent(in):: count

    integer(i8) :: i

    self%step = step
    log_debug(*) "Initializing coordinate '", trim(self%name), "' from complete range."

    ! set the bounds first
    select case(self%staggerId)
    case(0)
      self%bounds(1) = start
    case(1)
      self%bounds(1) = start - self%step / 2
    case(2)
      self%bounds(1) = start - self%step
    end select
    self%bounds(2) = self%bounds(1) + self%step * count

    ! set values and count properties
    call self%set_count(count)
    allocate(self%values(0_i8 : self%count))
    self%values(:) = self%bounds(1) + (/(i * self%step, i = 0_i8, self%count)/)

  end subroutine from_range

  subroutine from_other(self, other)
    !< initialize a Coordinate instance based on another Coordinate instance
    !< usually some properties are already set but some are unknown, these are taken from other Coordinate

    !> Coordinate type-bound procedure
    class(Coordinate), intent(inout) :: self
    !> other Coordinate to use remaining properties from
    class(Coordinate), intent(in) :: other

    integer(i4) :: i, oppositeBoundIndex
    real(dp) :: newBound

    ! okay so here we are not finalized and ...
    if (allocated(self%values)) then
      ! ... values is allocated, we have to set the remaining bound
      ! get remaining bound
      oppositeBoundIndex = self%get_bound_index()
      newBound = other%bounds(oppositeBoundIndex)
      ! set and check
      call self%insert_bound(newBound)
    else if (ne(self%step, nodata_dp) .or. self%count /= defaultCoordCount) then
      ! ... step or count is set so we need other step
      log_debug(*) "Initializing coordinate '", trim(self%name), "' from other coordinate '", trim(other%name), "'."
      if (self%count /= defaultCoordCount) then
        ! get the bounds from other
        self%step = (other%bounds(2) - other%bounds(1)) / self%count

      else if (ne(self%step, nodata_dp)) then
        if (eq(other%step, nodata_dp)) then
          log_error(*)  "The coordinate '", trim(self%name), "' cannot be initialized from other coordinate '", &
                  trim(other%name), "', as it does not have a constant step size. Other coordinate: "
          call other%get_stats()
          stop 1
        end if
        call self%set_count(floor(abs(other%count * other%step / self%step), kind=i8))
      end if
      allocate(self%values(0_i8: self%count))
      self%values = other%bounds(1) + (/(i * self%step, i = 0_i8, self%count)/)
      self%bounds(:) = [self%values(0), self%values(self%count)]
    else
      log_error(*)  "The coordinate '", trim(self%name), "' cannot be initialized, please set", &
              " a step size in the file ", filenameNamelistMprDefault
      stop 1
    end if

    log_debug(*) 'from_other: resulting values', self%values
    log_debug(*) 'from_other: resulting bounds', self%bounds
    log_debug(*) 'from_other: resulting step', self%step

  end subroutine from_other

  subroutine insert_bound(self, newBound)
    !< insert a missing bound in the existing temporary property of self%values
    !< this reallocates the values

    !> Coordinate type-bound procedure
    class(Coordinate), intent(inout) :: self
    !> the bound to be inserted
    real(dp) :: newBound

    real(dp), dimension(self%count) :: tempValues

    select case(self%staggerId)
      case(0)
        tempValues = self%values(0_i8 : self%count - 1_i8)
        deallocate(self%values)
        call self%set_values(tempValues, newBound)
      case(1)
        ! reverse the operations done at set_values
        ! this case should normally not be called, as the set_values routine issues an error
        ! if it is called without a bound and this routine should be called only when there was no bound available
        tempValues(2_i8:self%count) = self%values(1_i8:self%count-1_i8) - &
                (self%values(2_i8:self%count) - self%values(1_i8:self%count-1_i8)) / 2.0_dp
        tempValues(1) = self%values(1) - (self%values(2) - self%values(1)) / 2.0_dp
        deallocate(self%values)
        call self%set_values(tempValues, newBound)
        ! this warning is issued as vector values are manipulated if
        log_warn("(1X,A,A,F0.6,A,A,A)") "Bound (", newBound, ") for coordinate '", &
                trim(self%name), "' should not be re-inserted as it is center-staggered."
        ! stop 1
      case(2)
        tempValues = self%values(1_i8:self%count)
        deallocate(self%values)
        call self%set_values(tempValues, newBound)
    end select

  end subroutine insert_bound

  function get_bound_index(self) result(bound)
    !< returns the the index of the missing bound value needed for proper values initialization from self%bounds

    !> Coordinate type-bound procedure
    class(Coordinate), intent(in) :: self
    !> index of the missing bound value needed for proper values initialization from self%bounds
    integer(i4) :: bound

    select case(self%staggerId)
    case(0)
      bound = 2
    case(1:2)
      bound = 1
    end select

  end function get_bound_index

  function get_bound(self) result(bound)
    !< returns the missing bound value needed for proper values initialization

    !> Coordinate type-bound procedure
    class(Coordinate), intent(in) :: self
    !> the missing bound value needed for proper values initialization
    real(dp) :: bound

    select case(self%staggerId)
    case(0)
      bound = self%bounds(2)
    case(1:2)
      bound = self%bounds(1)
    end select

  end function get_bound

  subroutine set_bound(self, bound)
    !< sets the bound property as set during initialization

    !> Coordinate type-bound procedure
    class(Coordinate), intent(inout) :: self
    !> bound of values, upper bound of last cell, if start-staggered,
    !> lower bound of first cell if center or end-staggered
    real(dp), intent(in) :: bound

    select case(self%staggerId)
    case(0)
      self%bounds(2) = bound
    case(1:2)
      self%bounds(1) = bound
    end select

  end subroutine set_bound

  subroutine set_values(self, values, bound)
    !< sets the values property (ranges from 0 to self%count)
    !< index 0: lower bound of first cell
    !< index n: upper bound of nth cell

    !> Coordinate type-bound procedure
    class(Coordinate), intent(inout) :: self
    !> vector of values (can be either start-, center- or end-staggered)
    real(dp), dimension(:), intent(in) :: values
    !> optionally provide bound of values (start-staggered: bound is upper bound of last cell,
    !> center and end-staggered: bound is lower bound of first cell)
    !> if a step size is defined (regular increment), then bound is not needed
    real(dp), intent(in), optional :: bound

    integer(i8) :: i

    ! initialize from values
    call self%set_count(size(values, kind=i8))
    allocate(self%values(0_i8 : self%count))

    ! set the values depending on self%staggerId
    select case(self%staggerId)
    case(0)
      self%values(0_i8 : self%count - 1_i8) = values(:)
      ! start-bound values with bound determining the end boundary of the values
      if (present(bound)) then
        call self%check_bound(bound, values)
        self%values(self%count) = bound
      else if (ne(self%step, nodata_dp)) then
        ! use the step size
        self%values(self%count) = self%values(self%count - 1_i8) + self%step
      end if
    case(1)
      ! center-bound values with bound determining the start boundary of the values
      if (present(bound)) then
        call self%check_bound(bound, values)
        self%values(0) = bound
        ! center-bound values with different cell widths, need step-wise setting of cell sizes
        do i = 1_i8, self%count - 1_i8
          ! use the middle between the current and the following cell as the endpoint
          self%values(i) = values(i) + (values(i + 1_i8) - values(i)) / 2.0_dp
        end do
        if (self%count == 1_i8) then
          ! infer the end of the last cell from previous cell's width, values(i-1) does not exist in this case
          self%values(self%count) = values(i) + (values(i) - self%values(0)) / 2.0_dp
        else
          ! infer the end of the last cell from previous cell's width
          self%values(self%count) = values(i) + (values(i) - values(i - 1_i8)) / 2.0_dp
        end if
      else if (ne(self%step, nodata_dp)) then
        ! use the step size
        self%values(0_i8 : self%count - 1_i8) = values(:) - (self%step / 2.0_dp)
        self%values(self%count) = self%values(self%count - 1_i8) + self%step
      else
        log_error(*) 'Trying to set coordinate values with center-staggerd values and without a bound.', &
             'Unexpected configuration, contact developer!'
        call self%get_stats()
        stop 1
      end if
    case(2)
      ! end-bound values with bound determining the start boundary of the values
      self%values(1:) = values(:)
      if (present(bound)) then
        call self%check_bound(bound, values)
        self%values(0) = bound
      else if (ne(self%step, nodata_dp)) then
        ! use the step size
        self%values(0) = self%values(1) - self%step
      else if (self%count > 1_i8) then
        self%values(0) = self%values(1) - (self%values(2) - self%values(1))
      else
        log_error(*) 'Unexpected configuration, contact developer!'
        call self%get_stats()
        stop 1
      end if
    end select

    ! set bounds property
    self%bounds(1) = self%values(0)
    self%bounds(2) = self%values(self%count)
    log_debug(*) 'set_values: set values for Coordinate "', trim(self%name), '" with stagger "', &
            trim(self%staggerName), '" and length ', self%count, ' and bounds ', &
            self%bounds(1), ' and ', self%bounds(2), ' with flag isAscending ', self%is_ascending()


  end subroutine set_values

  subroutine check_bound(self, bound, values)
    !> Coordinate type-bound procedure
    class(Coordinate), intent(inout) :: self
    !> vector of values (can be either start-, center- or end-staggered)
    real(dp), dimension(:), intent(in) :: values
    !> optionally provide bound of values (start-staggered: bound is upper bound of last cell,
    !> center and end-staggered: bound is lower bound of first cell)
    !> if a step size is defined (regular increment), then bound is not needed
    real(dp), intent(in) :: bound

    integer(i8) :: isAscending
    character(512) :: errorString

    ! check for ascending values
    isAscending = 0_i8
    if (values(1) > values(size(values, kind=i8))) then
      isAscending = 1_i8
    end if
    ! prepare first part of possible error message
    errorString = "The provided values bound ("//compress(num2str(bound))//") for coordinate '"//trim(self%name)//"'"

    select case(self%staggerId + (3_i8 * isAscending))
    case(0)
      ! ascending and start-bound
      if (bound <= values(self%count)) then
        if (size(values, kind=i8) > 1_i8) then
          log_error("(1X,A,A,I0,A,A,A)") trim(errorString), " must be greater than the last value (", &
              values(self%count), ") if the values are ", trim(self%staggerName), "-staggered and in ascending order."
          stop 1
        end if
      end if
    case(1,2)
      ! ascending and center- or end-bound
      if (bound >= values(1)) then
        if (size(values, kind=i8) > 1_i8) then
          log_error("(1X,A,A,I0,A,A,A)") trim(errorString), " must be smaller than the first value (", &
              values(1), ") if the values are ", trim(self%staggerName), "-staggered and in ascending order."
          stop 1
        end if
      end if
    case(3)
      ! descending and start-bound
      if (bound >= values(self%count)) then
        if (size(values, kind=i8) > 1_i8) then
          log_error("(1X,A,A,I0,A,A,A)") trim(errorString), " must be smaller than the last value (", &
            values(self%count), ") if the values are ", trim(self%staggerName), "-staggered and in descending order."
          stop 1
        end if
      end if
    case(4,5)
      ! descending and center- or end-bound
      if (bound <= values(1)) then
        if (size(values, kind=i8) > 1_i8) then
          log_error("(1X,A,A,I0,A,A,A)") trim(errorString), " must be greater than the first value (", &
            values(1), ") if the values are ", trim(self%staggerName), "-staggered and in descending order."
          stop 1
        end if
      end if
    end select

  end subroutine check_bound

  subroutine set_step(self, values, step)
    !< sets the step property (difference between the cells)
    !< if irregular, step size equals nodata_dp

    !> Coordinate type-bound procedure
    class(Coordinate), intent(inout) :: self
    !> optionally provide a vector of values, this one is checked for whether the distances between
    !> two consecutive cells are equal within a tolerance (maxTolerance)
    real(dp), intent(in), dimension(:), optional :: values
    !> provide the step size directly
    real(dp), intent(in), optional :: step

    integer(i8) :: i, count

    if (present(step)) then
      self%step = step
      log_debug(*) 'set_step: set step for coordinate "', trim(self%name), '" based on step size passed: ', step
    else if (present(values)) then
      count = size(values, kind=i8)

      if (count >= 2_i8) then
        ! check if all values have similar length
        self%step = values(2) - values(1)
        do i = 2_i8, count
          ! if (ne((values(i) - values(i - 1)), self%step)) then
          ! TODO: how to handle this, what is numerical precision in netcdf dimensions?
          if (abs((values(i) - values(i - 1_i8)) - self%step) >= maxTolerance) then
            self%step = nodata_dp
            log_debug(*) 'set_step: set step for coordinate "', trim(self%name), '" to ', nodata_dp, &
                    'as there are differing step sizes (', values(i), ', ', values(i - 1), ') at index ', i
            exit
          end if
        end do
        log_debug(*) 'set_step: set step for coordinate "', trim(self%name), '" to ', self%step
      else
        log_debug(*) 'set_step: set step for coordinate "', trim(self%name), '" to ', nodata_dp, &
              'as the array has length 1 only'
        self%step = nodata_dp
      end if
    else
      log_error(*) 'Unexpected configuration, contact developer!'
      call self%get_stats()
      stop 1
    end if

  end subroutine set_step

  logical function is_ascending(self) result(isAscending)
    !< gets the isAscending property

    !> Coordinate type-bound procedure
    class(Coordinate), intent(in) :: self

    isAscending = defaultCoordIsAscending
    if (self%is_finalized()) then
      ! check for ascending values
      if (self%bounds(1) > self%bounds(2)) then
        isAscending = .false.
      end if
    else if (ne(self%step, nodata_dp)) then
      isAscending = self%step > 0.0_dp
    end if

  end function is_ascending

  subroutine set_count(self, count)
    !< sets the count property (number of cells, not the number of boundaries)

    !> Coordinate type-bound procedure
    class(Coordinate), intent(inout) :: self
    !> number of cells
    integer(i8), intent(in) :: count

    self%count = abs(count)
  end subroutine set_count

  subroutine set_2d_count(self)
    !< attempts to set the count property of 2d coordinates
    class(Coordinate), intent(inout) :: self

    integer(i4) :: iCoord
    integer(i4), dimension(:), allocatable :: coordIndex

    if (self%is_polygon()) then
      ! assume that last index of polygons is the number of polygons
      call self%set_count(size(self%centersCoord1, kind=i8))
    else if (self%is_2d()) then
      if (allocated(self%subDimSizes)) then
        deallocate(self%subDimSizes)
      end if
      allocate(self%subDimSizes(size(self%subDims)))
      ! check if the regular 2d coordinate variable is initialized
      do iCoord=1, size(self%subDims)
        coordIndex = get_index_in_coordinate(self%subDims(iCoord))
        ! should not have size 2
        self%subDimSizes(iCoord) = MPR_COORDINATES(coordIndex(1))%count
      end do
      call self%set_count(product(self%subDimSizes))
    end if

  end subroutine set_2d_count

  subroutine set_polygons_from_2d(self, centersOnlyArg)
    !< attempts to set the polygon properties out of regular 2d coordinates
    class(Coordinate), intent(inout) :: self
    logical, intent(in), optional :: centersOnlyArg

    logical :: centersOnly
    integer(i4) :: iCoord, iCorner
    integer(i8) :: iCell, iCell1, iCell2
    integer(i8), dimension(4), parameter :: coord1Addon = [0_i4, 1_i4, 1_i4, 0_i4]
    integer(i8), dimension(4), parameter :: coord2Addon = [0_i4, 0_i4, 1_i4, 1_i4]
    type(CoordinatePointer), dimension(:), allocatable :: sourceCoords

    ! set the pointers to the two source coordinates
    log_debug(*) 'set_polygons_from_2d: Calculating the cell corners and centers for coordinate ', trim(self%name)

    if (present(centersOnlyArg)) then
      centersOnly = centersOnlyArg
    else
      centersOnly = .false.
    end if
    ! check for rank
    if (self%rank /= 2_i4) then
      log_error(*) 'set_polygons_from_2d: Attempting to calculate polygons for a ', self%rank, &
      'D variable, only 2D supported.'
      stop 1
    end if
    self%corners =  2_i4 ** size(self%subDims)
    allocate(sourceCoords(self%rank))
    do iCoord=1, self%rank
      call sourceCoords(iCoord)%set_coordinate_pointer(self%subDims(iCoord))
    end do

    if (self%count == defaultCoordCount) then
      call self%set_2d_count()
    end if

    if (.not. centersOnly) then
      allocate(self%cornersCoord1(self%corners, self%count))
      allocate(self%cornersCoord2(self%corners, self%count))
    end if
    allocate(self%centersCoord1(self%count))
    allocate(self%centersCoord2(self%count))

    ! loop over each cell and set the midpoints and the corners
    do iCell=1_i8, self%count
      ! get the index of both source coordinates for the current target cell
      iCell1 = mod(iCell -1_i8, self%subDimSizes(1)) + 1_i8
      iCell2 = ((iCell -1_i8) / self%subDimSizes(1)) + 1_i8
      if (.not. centersOnly) then
        ! loop over all 4 corners and set the corners of the cell in a counterclockwise fashion (see dimxAddon)
        do iCorner=1, self%corners
          self%cornersCoord1(iCorner, iCell) = sourceCoords(1)%coord_p%values(iCell1-1_i8 + coord1Addon(iCorner))
          self%cornersCoord2(iCorner, iCell) = sourceCoords(2)%coord_p%values(iCell2-1_i8 + coord2Addon(iCorner))
        end do
      end if
      ! now set the center points, simply use mean of bounds for that
      self%centersCoord1(iCell) = (sourceCoords(1)%coord_p%values(iCell1-1_i8) + sourceCoords(1)%coord_p%values(iCell1)) / 2.0_dp
      self%centersCoord2(iCell) = (sourceCoords(2)%coord_p%values(iCell2-1_i8) + sourceCoords(2)%coord_p%values(iCell2)) / 2.0_dp
    end do


  end subroutine set_polygons_from_2d

  recursive function is_finalized(self)
    !< checks, whether the Coordinate type is fully initialized (vector of values is allocated)

    !> Coordinate type-bound procedure
    class(Coordinate), intent(in) :: self
    logical :: is_finalized
    integer(i4) :: iCoord
    integer(i4), dimension(:), allocatable :: coordIndex

    if (.not. self%is_initialized) then
      ! check if the polygon based coordinate is initialized
      is_finalized = .false.
      return
    else if (self%is_polygon()) then
      is_finalized = .true.
      ! check if the polygon based coordinate is initialized
      if (size(self%centersCoord1) == 0_i4 .or. size(self%centersCoord2) == 0_i4) then
        is_finalized = .false.
        return
      end if
    else if (self%is_2d()) then
      is_finalized = .true.
      ! check if the regular 2d coordinate variable is initialized
      do iCoord=1, size(self%subDims)
        coordIndex = get_index_in_coordinate(self%subDims(iCoord))
        ! should not have size 2
        if (.not. MPR_COORDINATES(coordIndex(1))%is_finalized()) then
          is_finalized = .false.
          return
        end if
      end do
    else
      ! check if the regular 1d coordinate variable is initialized
      is_finalized = allocated(self%values)
      if (is_finalized) then
        is_finalized = ne(self%values(0), nodata_dp) .and. ne(self%values(self%count), nodata_dp)
      end if
    end if

  end function is_finalized

  subroutine reset(self)
    !< restores all properties to its default, e.g. deallocating the dummyVector and all allocatable properties

    !> Coordinate type-bound procedure
    class(Coordinate), intent(inout) :: self

    if (self%is_initialized) then
      self%name = ''
      self%is_initialized = .false.
      if (allocated(self%values)) deallocate(self%values)
      if (allocated(self%subDims)) deallocate(self%subDims)
      if (allocated(self%subDimSizes)) deallocate(self%subDimSizes)
      if (allocated(self%cornersCoord1)) deallocate(self%cornersCoord1)
      if (allocated(self%cornersCoord2)) deallocate(self%cornersCoord2)
      if (allocated(self%centersCoord1)) deallocate(self%centersCoord1)
      if (allocated(self%centersCoord2)) deallocate(self%centersCoord2)
    end if

  end subroutine reset

  subroutine write_coordinate(self, nc, NcDim)
    class(Coordinate), intent(in) :: self
    class(NcDataset), intent(inout) :: nc
    type(NcDimension), intent(out) :: NcDim

    character(maxNameLength), dimension(size(self%attributes)) :: attributeNames
    character(maxStringLength), dimension(size(self%attributes)) :: attributeValuesChar
    logical, dimension(size(self%attributes)) :: maskIsChar
    integer(i4) :: iAtt

    maskIsChar = .false.
    ! select now all the character valued attributes
    do iAtt=1_i4, size(self%attributes)
      if (self%attributes(iAtt)%get_type() == 2_i4) then
        attributeNames(iAtt) = self%attributes(iAtt)%get_name()
        attributeValuesChar(iAtt) = self%attributes(iAtt)%get_value()
        maskIsChar(iAtt) = .true.
      end if
    end do
    if (.not. self%is_2d()) then
      NcDim = nc%setCoordinate(trim(self%name), int(self%count, kind=i4), self%values(:), self%staggerID, &
              pack(attributeNames, maskIsChar), pack(attributeValuesChar, maskIsChar))
    else if (self%is_polygon()) then
      NcDim = nc%setCoordinate(trim(self%name), int(self%count, kind=i4), &
              attribute_names=pack(attributeNames, maskIsChar),&
              attribute_values=pack(attributeValuesChar, maskIsChar), &
              centersDim1=self%centersCoord1, centersDim2=self%centersCoord2, &
              cornersDim1=self%cornersCoord1, cornersDim2=self%cornersCoord2, &
              subDimSizes=int(self%subDimSizes, kind=i4), units=self%unit)
    else
      call self%get_stats()
      log_error(*) "The 2-dimensional coordinate variable ", trim(self%name), " cannot be written to a netcdf file"
      stop 1
      !coord1 = get_coordinate(self%subDims(1))
      !coord2 = get_coordinate(self%subDims(2))
      !NcDim = nc%set2dCoordinateVariable(trim(self%name), self%count, self%values(:), self%staggerID, &
      !            attribute_names(maskIsChar), attribute_values(maskIsChar))
    end if

  end subroutine write_coordinate

  subroutine get_stats(self)
    !< prints some information on the Coordinate
    !< useful for debugging or printing before raising error messages

    !> Coordinate type-bound procedure
    class(Coordinate), intent(in) :: self

    integer(i4) :: iSubCoord

    log_info("(1X,A,A,A,A,I0,A,F0.6,A,F0.6,A,I0,A,A,A,L1)") 'Coordinate "', trim(self%name), '" and id: ',self%id, &
            ' with bounds: ', self%bounds(1), ' and ', self%bounds(2), ' and length: ', self%count, &
            ' and stagger ', trim(self%staggerName), ' and flag isAscending: ', self%is_ascending()
    if (eq(self%step, nodata_dp)) then
      log_info(*) 'The step size is irregular.'
    else
      log_info("(1X,A,A,F0.6)") 'The step size is ', self%step
    end if

    if (allocated(self%values)) then
      ! make some pretty output here
      if (self%count > 3_i8) then
        ! print values with index 1,2,...,last
        log_info("(1X,A,A,F0.6,A,F0.6,A,F0.6,A)") 'The (cell end) values are: [', self%values(1), &
                ', ', self%values(2), ', ..., ', self%values(self%count), ']'
      else if (self%count == 3_i8) then
        ! print values with index 1,2,last
        log_info("(1X,A,A,F0.6,A,F0.6,A,F0.6,A)") 'The (cell end) values are: [', self%values(1), &
                ', ', self%values(2), ', ', self%values(self%count), ']'
      else if (self%count == 2_i8) then
        ! print values with index 1,2
        log_info("(1X,A,A,F0.6,A,F0.6,A)") 'The (cell end) values are: [', self%values(1), ', ', self%values(2), ']'
      else
        log_info("(1X,A,A,F0.6,A)") 'The (cell end) value is: [', self%values(1), ']'
      end if
    else
      log_info(*) 'The values are not yet initialized.'
    end if

    if (self%is_polygon()) then
      log_info("(1X,A,A,I0,A,I0,A)") 'The Coordinate is polygon-based and contains ', &
              self%count, ' polygons with maximum', &
              self%corners, 'nodes.'
    else if (self%is_2d()) then
      log_info(*) 'The Coordinate is a regular 2d-coordinate variable with the subDims: '
      do iSubCoord=1, size(self%subDims)
          log_info(*) trim(self%subDims(iSubCoord))
      end do
    end if

  end subroutine get_stats

  function get_index_in_coordinate(coordinateName, coordinates, checkAliasesArg) result(ids)
    !< this function checks if the coordinateName is in a vector of Coordinates
    !< this vector can be provided as an optional dummy argument, else it falls back to MPR_COORDINATES
    !< optionally also checks for coordinate aliases of coordinateName (checkAliasesArg)
    !< in this case multiple indices can be returned (a coordinate can be in multiple aliases groups)

    !> coordinate name to find
    character(*), intent(in) :: coordinateName
    !> optional vector of CoordinatePointers to search for coordinateName
    type(CoordinatePointer), dimension(:), intent(in), optional:: coordinates
    !> flag whether to also check for aliases (usually .true. only if in combination with coordinates)
    logical, optional, intent(in) :: checkAliasesArg
    !> returns vector of indices (length 1 if checkAliasesArg == .false., else length >= 1)
    integer(i4), dimension(:), allocatable :: ids

    character(maxNameLength), dimension(:), allocatable :: aliases
    type(Coordinate), dimension(:), pointer, save :: globalVector => null()

    integer(i4) :: iCoordAlias, nCoords, iCoord, newIndex, requireAliases
    logical :: checkAliases
    integer(i4), dimension(:), allocatable :: coordinateIndices
    integer(i4) :: currentCoordIndex
    type(Coordinate), dimension(:), allocatable, target :: tempCoordinates

    checkAliases = .false.
    if (present(coordinates)) then
      allocate(coordinateIndices(size(coordinates)))
      do iCoord=1, size(coordinateIndices)
        ! get number of initialized coordinates
        nCoords = get_n_initialized(MPR_COORDINATES)
        ! get the correct index
        coordinateIndices(iCoord) = get_index_in_vector(coordinates(iCoord)%coord_p%name, MPR_COORDINATES)
        log_trace(*) 'get_index_in_coordinate: getting search coordinate ', &
                trim(coordinates(iCoord)%coord_p%name), ' from global list and got index ', &
                coordinateIndices(iCoord), ' (possible: ', nCoords, ')'
      end do
      ! this construct is needed as "Vector subscript array section in pointer assignment" is not allowed based
      ! on Fortran2008 standard
      tempCoordinates = MPR_COORDINATES(coordinateIndices)
      deallocate(coordinateIndices)
      globalVector => tempCoordinates
    else
      globalVector => MPR_COORDINATES
    end if

    ! get number of initialized coordinates
    nCoords = get_n_initialized(globalVector)
    ! allocate vector of size 0
    allocate(ids(0))
    ! get the index
    newIndex = get_index_in_vector(coordinateName, globalVector)

    if (present(checkAliasesArg)) then
      checkAliases = checkAliasesArg
    end if

    ! check for aliases and coordinateName not in vector?
    if (checkAliases .and. newIndex > nCoords) then
      ! check if coordinate has aliases attributed to it
      nCoords = get_n_initialized(MPR_COORDINATES)
      currentCoordIndex = get_index_in_vector(coordinateName, MPR_COORDINATES)
      if (currentCoordIndex > nCoords) then
        ! check if aliases exist
        call mpr_get_coord_alias(coordinateName, .false., aliases)
        requireAliases = 0_i4
      else
        ! the coordinate exists, now we pass additional argument to enable
        ! mpr_get_coord_alias providing sorted aliases depending on the order of MPR_COORDINATES(currentCoordIndex)%subDims
        call mpr_get_coord_alias(coordinateName, .false., aliases, MPR_COORDINATES(currentCoordIndex)%subDims)
        requireAliases = size(MPR_COORDINATES(currentCoordIndex)%subDims)
      end if
      log_trace(*) 'get_index_in_coordinate: could not find coordinate ', &
          trim(coordinateName), ', now checking for aliases: (', aliases, ')'

      ! check for aliases, now use the aliases as outer loop to implicitly consider a possible order of
      ! aliases as provided by mpr_get_coord_alias
      do iCoordAlias = 1, size(aliases)
        do iCoord = 1, size(globalVector)
          log_trace(*) 'get_index_in_coordinate: comparing ', trim(aliases(iCoordAlias)), ' and ', &
                  trim(globalVector(iCoord)%name)
          if (trim(aliases(iCoordAlias)) == trim(globalVector(iCoord)%name)) then
            ! append all the found indices
            call append(ids, iCoord)
          end if
        end do
      end do
      if (size(ids) < requireAliases) then
        log_error(*) 'Did not provide all ', requireAliases, &
                ' source coordinates for target coordinate ', trim(coordinateName), &
                '. I checked these aliases(', aliases, '). Check them in the section main - coordinate_aliases.'
        stop 1
      end if
      ! if no indexes were found, return foundIndex
      if (size(ids) == 0_i4) call append(ids, newIndex)
    else
      ! regularly append the foundIndex
      call append(ids, newIndex)
    end if

  end function get_index_in_coordinate

  subroutine add_coordinate(item, insertIndex)
    class(Coordinate), intent(in) :: item
    integer(i4), intent(in) :: insertIndex

    ! check if there is still room in the vector for one more coordinate
    if (insertIndex >= size(MPR_COORDINATES)) then
      log_error("(1X,A,A,I0,A,I0)")  "Coordinate cannot be added at index ('", insertIndex, &
              "'). Vector of coordinates only has length: ", size(MPR_COORDINATES)
      stop 1
    end if

    ! add item
    MPR_COORDINATES(insertIndex) = item

  end subroutine add_coordinate

  function create_upscaler_name_dimp(sourceCoords, targetCoords, upscaleOperatorNames) result(upscalerName)
    type(CoordinatePointer), dimension(:), intent(in) :: sourceCoords, targetCoords
    character(maxNameLength), dimension(:), intent(in), optional :: upscaleOperatorNames

    character(maxNameLength) :: upscalerName
    integer(i4) :: iCoord, iName

    upscalerName = ''
    ! get the source_dim names
    do iCoord = 1, size(sourceCoords)
      write(upscalerName, '(A,I0,A)') trim(upscalerName), sourceCoords(iCoord)%coord_p%id, '_'
    end do
    upscalerName = trim(upscalerName) // '_to__'

    ! get the target_dim names
    do iCoord = 1, size(targetCoords)
      write(upscalerName, '(A,I0,A)') trim(upscalerName), targetCoords(iCoord)%coord_p%id, '_'
    end do

    if (present(upscaleOperatorNames)) then
      upscalerName = trim(upscalerName) // '_with__'

      do iName=1, size(upscaleOperatorNames)
        upscalerName = trim(upscalerName) // trim(upscaleOperatorNames(iName)) // '_'
      end do

    end if

    ! remove trailing underscore
    upscalerName = upscalerName(:len_trim(upscalerName)-1)

    if (len_trim(upscalerName) >= maxNameLength) then
      log_error(*) 'attempting to create an upscaler name out of coordinate pointers,', &
              ' but coordinate names are too long: ', upscalerName
      stop 1
    end if

  end function create_upscaler_name_dimp

  function create_upscaler_name_dim(sourceCoords, targetCoords, upscaleOperatorNames) result(upscalerName)
    type(Coordinate), intent(in) :: sourceCoords, targetCoords
    character(maxNameLength), dimension(1), intent(in), optional :: upscaleOperatorNames
    character(maxNameLength) :: upscalerName

    ! get the source_dim names
    write(upscalerName, '(I0,A,I0)') sourceCoords%id, '__to__', targetCoords%id
    if (present(upscaleOperatorNames)) then
      upscalerName = trim(upscalerName) // '__with__' // trim(upscaleOperatorNames(1))
    end if

    if (len_trim(upscalerName) >= maxNameLength) then
      log_error(*) 'attempting to create an upscaler name out of coordinates,', &
              ' but coordinate names are too long: ', upscalerName
      stop 1
    end if

  end function create_upscaler_name_dim

  function create_upscaler_alias_dimp(sourceCoords, targetCoords) result(upscalerName)
    type(CoordinatePointer), dimension(:), intent(in) :: sourceCoords, targetCoords
    character(maxNameLength) :: upscalerName
    integer(i4) :: iCoord

    upscalerName = ''
    ! get the source_dim names
    do iCoord = 1, size(sourceCoords)
      upscalerName = trim(upscalerName) // trim(sourceCoords(iCoord)%coord_p%name) // '_'
    end do
    upscalerName = trim(upscalerName) // '_to__'

    ! get the target_dim names
    do iCoord = 1, size(targetCoords)
      upscalerName = trim(upscalerName) // trim(targetCoords(iCoord)%coord_p%name) // '_'
    end do
    ! remove trailing underscore
    upscalerName = upscalerName(:len_trim(upscalerName)-1)

    if (len_trim(upscalerName) >= maxNameLength) then
      log_error(*) 'attempting to create an upscaler name out of coordinate pointers,', &
              ' but coordinate names are too long: ', upscalerName
      stop 1
    end if

  end function create_upscaler_alias_dimp

  function create_upscaler_alias_dim(sourceCoords, targetCoords) result(upscalerName)
    type(Coordinate), intent(in) :: sourceCoords, targetCoords
    character(maxNameLength) :: upscalerName

    ! get the source_dim names
    upscalerName = trim(sourceCoords%name) // '__to__'  // trim(targetCoords%name)

    if (len_trim(upscalerName) >= maxNameLength) then
      log_error(*) 'attempting to create an upscaler name out of coordinates,', &
              ' but coordinate names are too long: ', upscalerName
      stop 1
    end if

  end function create_upscaler_alias_dim

  subroutine combine_coordinates(sourceCoords, id)
    !< combine multiple 1D coordinates to a 2D coordinate
    !> use a vector CoordinatePointers as input
    type(CoordinatePointer), dimension(:), intent(in) :: sourceCoords
    !> return the id in the global vector MPR_COORDINATES
    integer(i4), intent(out) :: id

    integer(i4), dimension(:), allocatable :: ids
    character(maxNameLength) :: newCoordName
    character(maxStringLength) :: unit
    integer(i4) :: nCoords, iCoord
    type(Coordinate) :: dim

    ! build new name
    if (size(sourceCoords) > 2_i4) then
      log_warn(*) 'combining more than 2 coordinates is still experimental and not tested thoroughly'
    end if
    newCoordName = '2D'
    do iCoord=1, size(sourceCoords)
      newCoordName = trim(newCoordName)//'_'//trim(sourceCoords(iCoord)%coord_p%name)
    end do
    ! get common unit
    unit = get_common_unit(sourceCoords(1)%coord_p, sourceCoords(2)%coord_p)
    log_trace(*) 'combining ', size(sourceCoords), ' coordinates to a new coordinate: ', trim(newCoordName)
    ! check if already created
    ids = get_index_in_coordinate(newCoordName)
    id = ids(1)
    nCoords = get_n_initialized(MPR_COORDINATES)
    if (id > nCoords) then
      ! create new dim
      dim = Coordinate(&
              newCoordName, &
              id, &
              'start', &
              subDims=[(sourceCoords(iCoord)%coord_p%name, iCoord=1, size(sourceCoords))], &
              unit=unit &
      )
      call add_coordinate(dim, id)
    end if

  end subroutine combine_coordinates

  subroutine create_integral_coordinate(sourceCoord, id)
    !< create an arbitrary coordinate of length 1
    !< which covers the whole range of the sourceCoord (used for broadcasting)
    !> use a Coordinate as input
    type(Coordinate), intent(in) :: sourceCoord
    !> return the id in the global vector MPR_COORDINATES
    integer(i4), intent(out) :: id

    integer(i4), dimension(:), allocatable :: ids
    character(maxNameLength) :: newCoordName
    integer(i4) :: nCoords
    type(Coordinate) :: dim

    newCoordName = '#'//trim(sourceCoord%name)//'_integral#'
    ! check if already created
    ids = get_index_in_coordinate(newCoordName)
    id = ids(1)
    nCoords = get_n_initialized(MPR_COORDINATES)
    if (id > nCoords) then
      ! create new dim
      dim = Coordinate(&
              newCoordName, &
              id, &
              'start', &
              values=[sourceCoord%bounds(1)], &
              bound=sourceCoord%bounds(2), &
              unit=sourceCoord%unit &
      )
      call add_coordinate(dim, id)
      ! also register in MPR_COORD_ALIAS
      call mpr_add_coord_alias(newCoordName, trim(sourceCoord%name))
    end if

  end subroutine create_integral_coordinate

  subroutine split_coordinate(dim, ids)
    !< splits a 2D coordinate into multiple 1D coordinates
    !> use a CoordinatePointer as input
    type(CoordinatePointer), intent(in) :: dim
    !> return the ids in the global vector of MPR_COORDINATES
    integer(i4), dimension(:), allocatable, intent(out) :: ids

    integer(i4), dimension(:), allocatable :: dimIds
    character(maxNameLength) :: dimName
    integer(i4) :: iCoord, id, nCoords
    !type(Coordinate) :: dim

    if (.not. dim%coord_p%is_2d() .or. dim%coord_p%is_polygon()) then
      log_error(*) 'split_coordinate: Splitting of coordinates only supported for 2D coordinate variables'
      stop 1
    end if
    allocate(ids(size(dim%coord_p%subDims)))
    do iCoord=1, size(dim%coord_p%subDims)
      dimName = dim%coord_p%subDims(iCoord)
      dimIds = get_index_in_coordinate(dimName)
      id = dimIds(1)
      nCoords = get_n_initialized(MPR_COORDINATES)
      if (id > nCoords) then
        log_error(*) 'split_coordinate: Cannot create new coordinate out of 2D coordinate variable'
        stop 1
        ! ! create new dim
        ! dim = Coordinate(&
        !         dimName, &
        !         id, &
        !         'start', &
        !         subDims=[(sourceCoords(iCoord)%coord_p%name, iCoord=1, size(sourceCoords))])
        ! call add_coordinate(dim, id)
      end if
      ids(iCoord) = id
    end do
    log_trace(*) 'split_coordinate: Split coordinate ', trim(dim%coord_p%name), ' into ', dim%coord_p%subDims
    ! check if already created

  end subroutine split_coordinate

  function get_common_unit(coord1, coord2) result(unit)
    type(Coordinate), intent(in) :: coord1, coord2
    character(maxNameLength) :: unit

    if (trim(coord1%unit) == trim(coord2%unit)) then
      unit = trim(coord1%unit)
    else
    log_warn(*) 'get_common_unit: Comparing unit of coordinate "', trim(coord1%name), '" and coordinate "', &
        trim(coord2%name), '" (', trim(coord1%unit), ', ', trim(coord2%unit), ') and setting them to :', trim(defaultAlias)
      unit = defaultAlias
    end if
  end function get_common_unit

  subroutine check_within_bounds(coord, bounds)
    !> Coordinate type-bound procedure
    type(Coordinate), intent(in) :: coord
    !> two bounds of other coordinate
    real(dp), dimension(2), intent(in) :: bounds

    if ((minval(bounds) - minval(coord%bounds)) > maxTolerance) then
      call coord%get_stats()
      log_error("(1X,A,A,A,A,F0.6,A,F0.6,A)") "For coordinate '", &
              trim(coord%name), "', its lowest value (", &
              minval(coord%bounds), ") is lower than the lower bound to check against (", &
              minval(bounds), ")."
      stop 1
    end if
    if ((maxval(coord%bounds) - maxval(bounds)) > maxTolerance) then
      call coord%get_stats()
      log_error("(1X,A,A,A,A,F0.6,A,F0.6,A)") "For coordinate '", &
              trim(coord%name), "', its highest value (", &
              maxval(coord%bounds), ") is higher than the upper bound to check against (", &
              maxval(bounds), ")."
      stop 1
    end if

  end subroutine check_within_bounds

  function has_attribute(coord, attName, attType, attLength) result(hasAttribute)
    !< \brief check if coordinate has attribute with name attName
    !< \details check if coordinate has attribute with name attName
    class(Coordinate), intent(in) :: coord
    character(*), intent(in) :: attName
    integer(i4), intent(out), optional :: attType
    integer(i4), intent(out), optional :: attLength
    logical :: hasAttribute

    character(maxNameLength) :: attNameCoord
    integer(i4) :: iAtt

    hasAttribute = .false.
    if (present(attType)) attType = nodata_i4
    if (present(attLength)) attLength = nodata_i4
    do iAtt = 1, size(coord%attributes)
      attNameCoord = coord%attributes(iAtt)%get_name()
      if (trim(attNameCoord) == trim(attName)) then
        hasAttribute = .true.
        if (present(attType)) attType = coord%attributes(iAtt)%get_type()
        if (present(attLength)) attLength = coord%attributes(iAtt)%get_length()
        return
      end if
    end do

  end function has_attribute

  subroutine get_attribute(coord, attName, attValue)
    !< \brief get coordinate's attribute value for given attribute name
    !< \details get coordinate's attribute value for given attribute name
    class(Coordinate), intent(in) :: coord
    character(*), intent(in) :: attName
    character(maxStringLength) :: attValue
    character(maxNameLength) :: attNameCoord
    integer(i4) :: attType

    integer(i4) :: iAtt

    attValue = ''
    do iAtt = 1, size(coord%attributes)
      attNameCoord = coord%attributes(iAtt)%get_name()
      attType = coord%attributes(iAtt)%get_type()
      if (trim(attNameCoord) == trim(attName) .and. attType == 2_i4) then
        attValue = coord%attributes(iAtt)%get_value()
        return
      end if
    end do

  end subroutine get_attribute

end module mo_mpr_coordinate

