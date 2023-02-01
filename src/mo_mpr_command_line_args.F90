!> \file mo_mpr_reset.f90
!> \brief reset all global lists of MPR types
!> \details reset all global lists of MPR types
!> \authors Robert Schweppe
!> \date Apr 2020

#include "flogging.h"
#ifndef MPRVERSION
#define MPRVERSION "0.0.0-dev0"
#endif
#ifndef MPRDATE
#define MPRDATE __DATE__
#endif
#define set_version(x) CHARACTER(len = *), PARAMETER :: version = x
#define set_date(x) CHARACTER(len = *), PARAMETER :: version_date = x

module mo_mpr_command_line_args

  use mo_kind, only: i4
  use mo_string_utils, only: splitString, replace_text
  use mo_mpr_constants, only: maxNameLength
  use flogging
  use mo_mpr_file, only: filenamenamelistmprdefault, filenamenamelistparamdefault

  implicit none

  set_version(MPRVERSION)
  !< Current mHM model version (will be set to \htmlinclude version.txt \latexinclude version.txt)

  set_date(MPRDATE)
  !< Time of current mHM model version release
  private

  public :: check_mpr_command_line_args

contains
  !> Check the command-line arguments to set the correct namelist paths
  subroutine check_mpr_command_line_args(filenameNamelistMpr, filenameNamelistParam, command)
    !> the paths to look for
    character(len=*), intent(inout) :: filenameNamelistMpr, filenameNamelistParam
    !> this arg is purely optional for testing purposes, replaces the command line
    character(len=*), intent(in), optional :: command

    integer(i4) :: i, status, count
    character(len=maxNameLength) :: arg
    ! character(len=maxNameLength) :: version_string
    character(len=maxNameLength), dimension(:), allocatable :: commands

    ! grab arguments from optional arg (testing purposes)
    if (present(command)) then
      commands = splitString(command, ' ')
      count = size(commands)
    else
    ! loop over each command line argument and perform some sanity checks
      count = command_argument_count()
      allocate(commands(count))
      do i=1, count
        call get_command_argument(i, arg, status=status)
        if (status == 0) then
          commands(i) = trim(arg)
        else
          log_error(*) 'check_mpr_command_line_args: could not properly retrieve command line argument #', &
          i, ': ', trim(arg)
          stop 1
        end if
      end do
    end if

    ! Loop over all command-line arguments to look for flags, then grab next argument
    do i=1, count
      select case (trim(commands(i)))
        case ("--parameter_file", "-p")
          filenameNamelistParam = trim(replace_text(replace_text(commands(i+1_i4), "'", ""), '"', ''))
        case ("--config_file", "-c")
          filenameNamelistMpr = trim(replace_text(replace_text(commands(i+1_i4), "'", ""), '"', ''))
        case ("--version")
          log_info(*) 'MPR version ', version, ' (', version_date, ')'
          stop 1
        case ("--help", "-h")
          log_info(*) 'MPR - Multiscale Parameter Regionalization (version ',  version, ') salutes you'
          log_info(*) 'available options:'
          log_info(*) '-c [--config_file] : path to configuration file, default: (', &
                  trim(filenamenamelistmprdefault),')'
          log_info(*) '-p [--parameter_file] : path to optional extra parameter configuration file, default: (', &
                  trim(filenamenamelistparamdefault),')'
          log_info(*) '-v[v[v]] : increase the verbosity of the output'
          log_info(*) '-q[q[q]] : decrease the verbosity of the output (-> quiet)'
          stop 1
      end select
    enddo

    deallocate(commands)

  end subroutine check_mpr_command_line_args

end module mo_mpr_command_line_args