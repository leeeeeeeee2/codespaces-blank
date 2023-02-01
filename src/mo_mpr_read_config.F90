!> \file mo_mpr_read_config.f90

!> \brief read mpr config

!> \details This module contains all mpr subroutines related to
!> reading the mpr configuration from file.

!> \authors Stephan Thober
!> \date Aug 2015
!         Modified, Robert Schweppe Dec 2017 - adapted for mpr

#include "flogging.h"
module mo_mpr_read_config

  use mo_kind, only : i4, dp, i8
  use flogging

  implicit none

  public :: mpr_read_config

contains

  subroutine mpr_read_config(file_namelist, unamelist, filenameNamelistParam, uNamelistParam, parameterValues, &
                  parameterNames)

    use mo_mpr_global_variables, only : &
            OUT_FILENAME, CHECK_FOR_NODATAVALUE, WRITE_WEIGHTS
    use mo_mpr_constants, only : &
            maxNoDataArrays, &
            maxNoInputFieldsPerDA, &
            maxNoParameterPerDA, &
            maxNoCoords, &
            maxNoCoordsPerDA, &
            maxNoCoordAliases, &
            maxCoordLen, &
            maxNoCoordUpscalers, &
            maxNoAttributes, &
            defaultAlias, &
            defaultCoordStagger, &
            defaultCoordCount, &
            defaultDataArrayLimits, &
            defaultDataArrayToFile, &
            defaultReadWeights, &
            defaultWriteWeights, &
            defaultSubcellIdsFieldName, &
            defaultWeightsFieldName, &
            defaultNSubcellsFieldName, &
            maxStringLength, maxNameLength, &
            defaultCheckforNoDataValue
    use mo_mpr_coordinate, only : Coordinate, add_coordinate, MPR_COORDINATES
    use mo_nml, only : open_nml, close_nml, position_nml
    use mo_constants, only : nodata_dp, nodata_i4
    use mo_mpr_util_types, only : MPR_COORD_ALIAS, CoordAlias
    use mo_utils, only : ne
    use mo_mpr_utils, only : get_index_in_vector
    use mo_mpr_data_array, only : add_data_array, DataArray, MPR_DATA_ARRAYS
    use mo_mpr_data_array_upscale, only : CoordUpscaler, add_coord_upscaler, MPR_COORD_UPSCALERS, &
            get_index_in_coord_upscaler
    use mo_mpr_parameters, only : MPR_PARAMETERS, mpr_add_parameters

    implicit none

    character(*), intent(in) :: file_namelist
    integer, intent(in) :: unamelist
    character(*), intent(in), optional :: filenameNamelistParam
    integer, intent(in), optional :: uNamelistParam
    real(dp), dimension(:), intent(in), optional :: parameterValues
    character(maxNameLength), dimension(:), intent(in), optional :: parameterNames

    ! variables for DataArray
    character(maxStringLength), dimension(maxNoDataArrays) :: from_file
    character(maxStringLength), dimension(maxNoDataArrays) :: transfer_func
    character(maxNameLength), dimension(maxNoDataArrays) :: name
    character(maxNameLength), dimension(maxNoDataArrays) :: transfer_func_label
    character(maxNameLength), dimension(maxNoInputFieldsPerDA, maxNoDataArrays) :: from_data_arrays
    real(dp), dimension(maxNoParameterPerDA, maxNoDataArrays) :: from_parameter_values
    character(maxNameLength), dimension(maxNoParameterPerDA, maxNoDataArrays) :: from_parameter_names
    character(maxNameLength), dimension(maxNoCoordsPerDA, maxNoDataArrays) :: target_coord_names
    character(maxNameLength), dimension(maxNoCoordsPerDA, maxNoDataArrays) :: upscale_ops
    real(dp), dimension(2, maxNoDataArrays) :: limits
    logical, dimension(maxNoDataArrays) :: to_file

    ! variables for Coordinate
    character(maxNameLength), dimension(maxNoCoordAliases, maxNoCoords) :: coordinate_aliases
    character(maxNameLength), dimension(maxNoCoords) :: coord_name
    character(6), dimension(maxNoCoords) :: coord_stagger
    character(maxStringLength), dimension(maxNoCoords) :: coord_from_file
    character(maxStringLength), dimension(maxNoCoords) :: coord_unit
    real(dp), dimension(maxNoCoords) :: coord_from_range_start
    integer(i8), dimension(maxNoCoords) :: coord_from_range_count
    real(dp), dimension(maxNoCoords) :: coord_from_range_step
    real(dp), dimension(maxNoCoords) :: coord_from_values_bound
    real(dp), dimension(maxCoordLen, maxNoCoords) :: coord_from_values
    character(maxNameLength), dimension(maxNoAttributes, maxNoCoords) :: coord_attribute_names
    character(maxStringLength), dimension(maxNoAttributes, maxNoCoords) :: coord_attribute_values
    character(maxStringLength), dimension(maxNoCoords) :: coord_proj_string
    character(maxNameLength), dimension(maxNoCoordAliases, maxNoCoords) :: coord_sub_dims

    ! variables for CoordUpscler
    character(maxNameLength), dimension(maxNoCoordUpscalers) :: upscaler_name
    character(maxStringLength), dimension(maxNoCoordUpscalers) :: from_weight_file
    character(maxNameLength), dimension(maxNoCoordUpscalers) :: subcell_ids_field_name, &
            weights_field_name, n_subcells_field_name

    integer(i4) :: i, j, nCoordAlias, id, k, item_counter
    type(Coordinate) :: coord_
    type(DataArray) :: dataarray_
    type(CoordUpscaler) :: coord_upscaler_
    logical :: read_weights
    logical, dimension(:), allocatable :: mask
    character(maxNameLength), dimension(:), allocatable :: tempTargetCoordNames
    character(maxNameLength), dimension(:), allocatable :: tempUpscaleOps
    character(maxNameLength), dimension(:), allocatable :: TempCoordsubDims
    logical :: doesFileExist

    namelist/main/ OUT_FILENAME, coordinate_aliases, WRITE_WEIGHTS, read_weights, CHECK_FOR_NODATAVALUE
    namelist/data_arrays/ name, from_file, from_data_arrays, transfer_func, transfer_func_label, target_coord_names, &
            upscale_ops, limits, to_file, from_parameter_values, from_parameter_names
    namelist/coordinates/ coord_name, coord_stagger, coord_from_file, coord_from_values, &
            coord_from_range_start, coord_from_range_step, coord_from_range_count, coord_from_values_bound, &
            coord_attribute_names, coord_attribute_values, &
            coord_proj_string, coord_sub_dims, coord_unit
    namelist/upscalers/ upscaler_name, from_weight_file, &
            subcell_ids_field_name, weights_field_name, n_subcells_field_name

    !===============================================================
    !  Initialization
    !===============================================================
    ! for dim aliases
    coordinate_aliases = defaultAlias
    ! for upscalers
    read_weights = defaultReadWeights
    WRITE_WEIGHTS = defaultWriteWeights
    upscaler_name = defaultAlias
    from_weight_file = defaultAlias
    subcell_ids_field_name = defaultSubcellIdsFieldName
    weights_field_name = defaultWeightsFieldName
    n_subcells_field_name = defaultNSubcellsFieldName
    ! for Coordinate
    coord_name = defaultAlias
    coord_stagger = defaultCoordStagger
    coord_from_file = defaultAlias
    coord_attribute_names = defaultAlias
    coord_attribute_values = defaultAlias
    coord_proj_string = defaultAlias
    coord_sub_dims = defaultAlias
    coord_from_range_start = nodata_dp
    coord_from_range_count = defaultCoordCount
    coord_from_range_step = nodata_dp
    coord_from_values_bound = nodata_dp
    coord_from_values = nodata_dp
    coord_unit = defaultAlias
    ! for DataArray
    name = defaultAlias
    upscale_ops = defaultAlias
    target_coord_names = defaultAlias
    from_data_arrays = defaultAlias
    from_parameter_names = defaultAlias
    from_parameter_values = nodata_dp
    from_file = defaultAlias

    do i = 1, maxNoDataArrays
      limits(:, i) = defaultDataArrayLimits(:)
    end do

    transfer_func = defaultAlias
    transfer_func_label = defaultAlias
    to_file = defaultDataArrayToFile
    CHECK_FOR_NODATAVALUE = defaultCheckforNoDataValue

    !===============================================================
    !  Read the (variable) parameters
    !===============================================================
    if (present(filenameNamelistParam) .and. present(uNamelistParam)) then
      ! check if file exists, if not skip those parameters
      inquire(FILE=filenameNamelistParam, EXIST=doesFileExist)
      if (doesFileExist) then
        call open_nml(filenameNamelistParam, uNamelistParam, quiet = .true.)
        call read_parameters('parameters', uNamelistParam)
        call close_nml(uNamelistParam)
      end if
    end if
    if (present(parameterValues) .and. present(parameterNames)) then
      call mpr_add_parameters(parameterNames(:), parameterValues(:))
    end if

    !===============================================================
    !  Read namelist specifying the model configuration
    !===============================================================
    call open_nml(file_namelist, unamelist, quiet = .true.)
    call position_nml('main', unamelist)
    read(unamelist, nml = main)

    !===============================================================
    !  Init the CoordAlias groups
    !===============================================================
    ! infer the number of needed CoordAlias groups
    nCoordAlias = 0_i4
    do i = 1, maxNoCoords
      ! use only initialized dim aliases
      mask = [(trim(coordinate_aliases(j, i)) /= trim(defaultAlias), j = 1, maxNoCoordAliases)]
      if (any(mask)) then
        nCoordAlias = nCoordAlias + 1_i4
      end if
    end do

    ! allocate the needed groups ...
    allocate(MPR_COORD_ALIAS(nCoordAlias))
    k = 1_i4
    ! ... and iteratively set them
    do i = 1, maxNoCoords
      ! use only initialized dim aliases
      mask = [(trim(coordinate_aliases(j, i)) /= trim(defaultAlias), j = 1, maxNoCoordAliases)]
      if (any(mask)) then
        MPR_COORD_ALIAS(k) = CoordAlias(pack(coordinate_aliases(:, i), mask))
        k = k + 1_i4
      end if
    end do
    deallocate(mask)

    !===============================================================
    !  Init the Coordinates
    !===============================================================
    call position_nml('coordinates', unamelist)
    read(unamelist, nml = coordinates)

    log_info(*) 'Initializing coordinates from file ', trim(file_namelist)
    item_counter = 0_i4
    do i = 1, maxNoCoords
      ! stop if no name is initialized
      if (trim(coord_name(i)) == trim(defaultAlias)) then
        cycle
      end if
      log_debug(*) 'Initializing coordinates: ', trim(coord_name(i))
      ! first get the index at which it should be inserted in global vector
      id = get_index_in_vector(coord_name(i), MPR_COORDINATES)

      TempCoordsubDims = pack(coord_sub_dims(:, i), &
                      ! use only initialized vector entries
                      [(trim(coord_sub_dims(j, i)) /= trim(defaultAlias), j = 1, maxNoCoordAliases)])
      ! ... init the derived type ...
      coord_ = Coordinate(&
              name=coord_name(i), &
              id=id, &
              stagger=coord_stagger(i), &
              fileName=coord_from_file(i), &
              values=pack(coord_from_values(:, i), &
                      ! use only initialized vector entries
                      ne(coord_from_values(:, i), nodata_dp)), &
              bound=coord_from_values_bound(i), &
              start=coord_from_range_start(i), &
              step=coord_from_range_step(i), &
              count_=coord_from_range_count(i), &
              attributeNames=pack(coord_attribute_names(:, i), &
                      ! use only initialized vector entries
                      [(trim(coord_attribute_names(j, i)) /= trim(defaultAlias), j = 1, maxNoAttributes)]), &
              attributeValuesChar=pack(coord_attribute_values(:, i), &
                      ! use only initialized vector entries
                      [(trim(coord_attribute_values(j, i)) /= trim(defaultAlias), j = 1, maxNoAttributes)]), &
              projString=coord_proj_string(i), &
              ! set the coord_sub_dims in inverse order (C vs. F ordering)
              subDims=[(TempCoordsubDims(j), j=size(TempCoordsubDims), 1, -1)], &
              unit=coord_unit(i) &
              )
      ! ... and add it to global vector
      call add_coordinate(coord_, id)
      item_counter = item_counter + 1_i4
    end do

    if (item_counter == 0_i4) then
      log_warn(*) 'You did not supply any Coordinates in the namelist.'
    else
      item_counter = 0_i4
      deallocate(TempCoordsubDims)
    end if

    !===============================================================
    ! Init the fixed parameters
    !===============================================================
    log_info(*) 'Initializing parameters from file ', trim(file_namelist)
    call read_parameters('parameters', unamelist)

    !===============================================================
    !  Init the CoordUpscalers
    !===============================================================
    if (read_weights) then
      call position_nml('upscalers', unamelist)
      read(unamelist, nml = upscalers)

      log_info(*) 'Initializing upscalers from file ', trim(file_namelist)
      do i = 1, maxNoCoordUpscalers
        ! stop if no name is initialized
        if (trim(upscaler_name(i)) == trim(defaultAlias)) then
          exit
        end if
        log_debug(*) 'Initializing upscaler: ', trim(upscaler_name(i))
        ! first get the index at which it should be inserted in global vector
        id = get_index_in_coord_upscaler(upscaler_name(i), MPR_COORD_UPSCALERS)
        ! ... init the derived type ...
        coord_upscaler_ = CoordUpscaler(&
                name = upscaler_name(i), &
                id = id, &
                alias = upscaler_name(i), &
                fromWeightFile = from_weight_file(i), &
                subcellIdsFieldName = subcell_ids_field_name(i), &
                weightsFieldName = weights_field_name(i), &
                nSubcellsFieldName = n_subcells_field_name(i) &
                )
        ! ... and add it to global vector
        call add_coord_upscaler(coord_upscaler_, id)
      end do
    end if

    !===============================================================
    !  Init the DataArrays
    !===============================================================
    call position_nml('data_Arrays', unamelist)
    read(unamelist, nml = data_arrays)

    log_info(*) 'Initializing data arrays from file ', trim(file_namelist)
    do i = 1, maxNoDataArrays
      ! stop if no name is initialized
      if (trim(name(i)) == trim(defaultAlias)) then
        cycle
      end if
      log_debug(*) 'Initializing data array: ', trim(name(i))
      ! first get the index at which it should be inserted in global vector
      id = get_index_in_vector(name(i), MPR_DATA_ARRAYS)

      ! check for parameter_names and replace them by the values
      do j = 1, maxNoParameterPerDA
        if (trim(from_parameter_names(j, i)) == trim(defaultAlias)) then
          exit
        end if
        if (ne(from_parameter_values(j, i), nodata_dp)) then
          log_error("(1X,A,A,A,A,I0)") &
                  "you supplied both 'from_parameter_names' and 'from_parameters_values' for the data array '", &
                  trim(name(i)), "' at index ", i
          stop 1
        end if
        if (MPR_PARAMETERS%does_contain(from_parameter_names(j, i))) then
          from_parameter_values(j, i) = MPR_PARAMETERS%get_parameter(from_parameter_names(j, i))
        else
          log_error("(1X,A,A,A,A,A,A,I0,A)") "the requested parameter name '", trim(from_parameter_names(j, i)), &
                  "' for data array '", trim(name(i)), "' at index ", i, " cannot be found."
          stop 1
        end if
      end do

      tempTargetCoordNames = pack(target_coord_names(:, i), &
                      ! use only initialized dim name
                      [(trim(target_coord_names(j, i)) /= trim(defaultAlias), j = 1, maxNoCoordsPerDA)])
      tempUpscaleOps = pack(upscale_ops(:, i), &
                      ! use only initialized upscale_op entries
                      [(trim(upscale_ops(j, i)) /= trim(defaultAlias), j = 1, maxNoCoordsPerDA)])

      ! ... init the derived type ...
      dataarray_ = DataArray(&
              name(i), &
              id, &
              pack(from_data_arrays(:, i), &
                      ! use only initialized field name
                      [(trim(from_data_arrays(j, i)) /= trim(defaultAlias), j = 1, maxNoInputFieldsPerDA)]), &
              from_file(i), &
              transfer_func(i), &
              transfer_func_label(i), &
              pack(from_parameter_values(:, i), &
                      ! use only initialized parameter name
                      [(ne(from_parameter_values(j, i), nodata_dp), j = 1, maxNoParameterPerDA)]), &
              ! set the target_coord_names in inverse order (C vs. F ordering)
              [(tempTargetCoordNames(j), j=size(tempTargetCoordNames), 1, -1)], &
              ! set the upscale_ops in inverse order (C vs. F ordering)
              [(tempUpscaleOps(j), j=size(tempUpscaleOps), 1, -1)], &
              limits(:, i), &
              to_file(i) &
              )
      ! ... and add it to global vector
      call add_data_array(dataarray_, id)
      item_counter = item_counter + 1_i4
    end do

    call close_nml(unamelist)

    if (item_counter == 0_i4) then
      log_error(*) 'You did not supply any DataArrays in the namelist.'
      stop 1
    else
      item_counter = 0_i4
      deallocate(tempUpscaleOps, tempTargetCoordNames)
    end if

  end subroutine mpr_read_config

  subroutine read_parameters(section_name, unamelist)
    use mo_mpr_parameters, only : &
            mpr_add_parameters
    use mo_mpr_constants, only : &
            maxNoParameters, &
            defaultAlias, maxNameLength
    use mo_nml, only : position_nml
    use mo_constants, only : nodata_dp

    character(*), intent(in) :: section_name
    integer, intent(in) :: unamelist

    character(maxNameLength), dimension(maxNoParameters) :: parameter_names
    real(dp), dimension(maxNoParameters) :: parameter_values

    !===============================================================
    ! Init the parameters
    !===============================================================
    namelist/parameters/ parameter_names, parameter_values

    ! for Parameters
    parameter_names = defaultAlias
    parameter_values = nodata_dp

    call position_nml(section_name, unamelist)
    read(unamelist, nml = parameters)

    call mpr_add_parameters(parameter_names(:), parameter_values(:))

  end subroutine read_parameters
end module mo_mpr_read_config
