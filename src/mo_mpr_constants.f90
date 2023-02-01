!> \file mo_mpr_constants.f90

!> \brief Provides MPR specific constants

!> \details Provides MPR specific constants such as flood plain elevation.

!> \author Matthias Cuntz
!> \date Nov 2011

MODULE mo_mpr_constants

  USE mo_kind, ONLY : i4, dp, i8
  use mo_constants, only: nodata_dp

  IMPLICIT NONE

  PRIVATE

  ! ------------------------------------------------------------------
  ! maximum values
  ! ------------------------------------------------------------------
  integer(i4), parameter, public :: maxStringLength = 2048_i4   ! maximum string length, e.g. for transfer functions
  integer(i4), parameter, public :: maxNameLength = 256_i4   ! maximum string length, e.g. for coordinate names etc.
  integer(i4), parameter, public :: maxNoDataArrays = 1000_i4     ! maximum number of data arrays
  integer(i4), parameter, public :: maxNoUpscalers = 50_i4     ! maximum number of upscalers
  integer(i4), parameter, public :: maxNoCoordUpscalers = 50_i4     ! maximum number of dim upscalers
  integer(i4), parameter, public :: maxNoParameters = 10000_i4     ! maximum number of constants
  integer(i4), parameter, public :: maxNoInputFieldsPerDA = 10_i4     ! maximum number of input fields per DA
  integer(i4), parameter, public :: maxNoParameterPerDA = 30_i4     ! maximum number of parameters per DA
  integer(i4), parameter, public :: maxNoCoordsPerDA = 5_i4     ! maximum number of allowed coordinates per DA
  integer(i4), parameter, public :: maxNoCoords = 30_i4     ! maximum number of allowed coordinates
  integer(i4), parameter, public :: maxNoCoordAliases = 10_i4     ! maximum number of coordinate aliases
  integer(i4), parameter, public :: maxCoordLen = 1000_i4     ! maximum length of coordinate
  integer(i4), parameter, public :: maxNoAttributes = 30_i4     ! maximum number of variable attributes
  integer(i4), parameter, public :: maxNoAttributeValues = 12_i4     ! maximum number of attribute values
  integer(i4), parameter, public :: maxNoTemporalData = 100_i4     ! maximum length of temporal data
  real(dp), parameter, public :: maxTolerance = 1.0e-9_dp     ! maximum tolerance for some equality comparisons

  ! ------------------------------------------------------------------
  ! default settings (for mpr.nml)
  ! ------------------------------------------------------------------
  character(6), parameter, public :: noname = '<NONE>' ! used in mo_mpr_reorder_data_array
  character(maxNameLength), parameter, public :: defaultAlias = "not_initialized"
  character(6), parameter, public :: defaultCoordStagger = "center"
  character(7), parameter, public :: defaultCoordUnits = "degrees"
  integer(i8), parameter, public :: defaultCoordCount = 0_i8
  logical, parameter, public :: defaultCoordIsAscending = .true.
  real(dp), parameter, public, dimension(2) :: defaultDataArrayLimits = nodata_dp
  logical, parameter, public :: defaultDataArrayToFile = .true.
  logical, parameter, public :: defaultReadWeights = .false.
  logical, parameter, public :: defaultWriteWeights = .false.
  logical, parameter, public :: defaultCheckforNoDataValue = .false.
  character(*), parameter, public :: defaultSubcellIdsFieldName = "subcell_ids"
  character(*), parameter, public :: defaultWeightsFieldName = "weights"
  character(*), parameter, public :: defaultNSubcellsFieldName = "n_subcells"
  character(12), parameter, public :: upperBoundName = ".upper_bound"
  character(12), parameter, public :: lowerBoundName = ".lower_bound"

  ! ------------------------------------------------------------------
  ! SCRIP conventions
  ! ------------------------------------------------------------------
  character(*), parameter, public :: defaultMapMethod = defaultAlias
  character(*), parameter, public :: scripSrcAddressName = 'src_address'
  character(*), parameter, public :: scripDstAddressName = 'dst_address'
  character(*), parameter, public :: scripWeightsName = 'remap_matrix'
  character(*), parameter, public :: scripWeightsCoord = 'num_wgts'
  character(*), parameter, public :: scripLinksCoord = 'num_links'

END MODULE mo_mpr_constants
