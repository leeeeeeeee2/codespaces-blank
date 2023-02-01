!> \file mo_mpr_util_types.f90

!> \brief MPR related types

!> \details A wrapper around the NetCDF Fortran 90 interface.
!
!> \authors Robert Schweppe
!> \date Jan 2018

#include "flogging.h"
module mo_mpr_util_types

  use mo_kind, only: i4
  use flogging
  use mo_mpr_constants, only : maxNameLength

  implicit none

  private

  public :: CoordAlias, mpr_get_coord_alias, mpr_add_coord_alias, MprBaseType, MPR_COORD_ALIAS

  ! --------------------------------------------------------------------------------------
  type :: CoordAlias

    character(maxNameLength), dimension(:), allocatable :: names

  contains
    procedure, public :: does_contain
    procedure, public :: get_alias
    procedure, public :: add_alias

  end type CoordAlias

  interface CoordAlias
    procedure newCoordAlias
  end interface CoordAlias

  type(CoordAlias), dimension(:), allocatable :: MPR_COORD_ALIAS ! all coord_alias

  type, abstract :: MprBaseType
    character(maxNameLength) :: name
    integer(i4) :: id
    logical :: is_initialized
  contains
    procedure(is_finalized_abstract), deferred :: is_finalized
    procedure(reset_abstract), deferred :: reset
  end type MprBaseType

  abstract interface
    function is_finalized_abstract(self) result(is_finalized)
      import MprBaseType
      class(MprBaseType), intent(in) :: self
      logical :: is_finalized
    end function is_finalized_abstract
  end interface

  abstract interface
    subroutine reset_abstract(self)
      import MprBaseType
      class(MprBaseType), intent(inout) :: self
    end subroutine reset_abstract
  end interface

  ! --------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------
contains

  function newCoordAlias(aliases)
    type(CoordAlias) :: newCoordAlias
    character(*), dimension(:), intent(in) :: aliases

    allocate(newCoordAlias%names(size(aliases)))
    newCoordAlias%names = aliases

  end function newCoordAlias

  function does_contain(self, name)
    class(CoordAlias), intent(in) :: self
    character(*), intent(in) :: name
    logical :: does_contain
    integer(i4) :: iName

    does_contain = any([(trim(self%names(iName)) == trim(name), iName=1, size(self%names))])
  end function does_contain

  subroutine get_alias(self, aliases, name)
    class(CoordAlias), intent(in) :: self
    character(*), dimension(:), allocatable, intent(inout) :: aliases
    character(*), intent(in), optional :: name

    if (present(name)) then
      aliases = pack(self%names, self%names /= trim(name))
    else
      aliases = self%names
    end if
  end subroutine get_alias

  subroutine add_alias(self, name)
    class(CoordAlias), intent(inout) :: self
    character(*), intent(in) :: name

    character(maxNameLength), dimension(:), allocatable :: tempNames

    if (.not. self%does_contain(name)) then
      allocate(tempNames(size(self%names) + 1))
      tempNames(size(self%names) + 1) = trim(name)
      tempNames(1:size(self%names)) = self%names
      call move_alloc(tempNames, self%names)
    end if

  end subroutine add_alias

  subroutine mpr_get_coord_alias(name, include_name, aliases, sortedAliasesFromGroup)
    !< get all the aliases of a coordinate, a coordinate can have multiple aliases and all these are collected
    character(*), intent(in) :: name
    logical, intent(in) :: include_name
    character(*), dimension(:), allocatable, intent(out) :: aliases
    character(*), dimension(:), allocatable, intent(in), optional :: sortedAliasesFromGroup

    character(maxNameLength), dimension(:), allocatable:: tempAliases
    character(maxNameLength), dimension(:), allocatable:: newAliases
    integer(i4) :: iCoordAlias, iAlias
    logical :: is_contained
    type(CoordAlias), dimension(:), allocatable :: aliasSelection

    if (present(sortedAliasesFromGroup)) then
      if (size(sortedAliasesFromGroup) > 0_i4) then
        ! select only the groups where sortedAliasesFromGroup are contained
        allocate(aliasSelection(size(sortedAliasesFromGroup)))
        do iAlias=1, size(sortedAliasesFromGroup)
          do iCoordAlias = 1, size(MPR_COORD_ALIAS)
            is_contained = MPR_COORD_ALIAS(iCoordAlias)%does_contain(sortedAliasesFromGroup(iAlias))
            if (is_contained) then
              aliasSelection(iAlias) = MPR_COORD_ALIAS(iCoordAlias)
            end if
          end do
        end do
        log_trace(*) 'mpr_get_coord_alias: provided, ', size(sortedAliasesFromGroup), &
              ' aliasesGroups to sort seach for coordinate : ', trim(name)
      else
        aliasSelection = MPR_COORD_ALIAS
      end if
    else
      ! select all aliasGroups
      aliasSelection = MPR_COORD_ALIAS
    end if

    allocate(aliases(0))
    do iCoordAlias = 1, size(aliasSelection)
      is_contained = aliasSelection(iCoordAlias)%does_contain(name)
      if (is_contained) then
        if (include_name) then
          call aliasSelection(iCoordAlias)%get_alias(newAliases)
        else
          call aliasSelection(iCoordAlias)%get_alias(newAliases, name)
        end if
        ! append the newAliases at the end of the existing aliases
        allocate(tempAliases(size(aliases)+size(newAliases)))
        tempAliases(1:size(aliases)) = aliases(:)
        tempAliases(size(aliases)+1:) = newAliases(:)
        call move_alloc(tempAliases, aliases)
      end if
    end do

  end subroutine mpr_get_coord_alias

  subroutine mpr_add_coord_alias(new_name, alias)
  !< add a new coordinate name (new_name) to an existing MPR_COORD_ALIAS group, if it matches alias (functions as key)
    character(*), intent(in) :: new_name
    character(*), intent(in) :: alias

    integer(i4) :: iCoordAlias
    logical :: is_contained

    do iCoordAlias = 1, size(MPR_COORD_ALIAS)
      is_contained = MPR_COORD_ALIAS(iCoordAlias)%does_contain(alias)
      if (is_contained) then
        call MPR_COORD_ALIAS(iCoordAlias)%add_alias(new_name)
        exit
      end if
    end do

  end subroutine mpr_add_coord_alias

end module mo_mpr_util_types

