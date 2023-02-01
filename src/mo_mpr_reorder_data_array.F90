!> \file mo_mpr_reorder_data_array.f90

!> \brief Reorder MPR_data_array

!> \details Optimize order in which data arrays are executed to minimize the MPR storage requirement
!>          and make order of execution independent of order provided by the user in the mpr.nml. The
!>          reordering assures that all required input fields are available at the time that they are needed.
!>          The algorithm to select a node works like this:
!>          1) algorithm always picks element with depth one.
!>          2) if there are multiple, it selects the one where most previously calculated data arrays can be deleted
!>          3) if nothing can be deleted, it picks the one with the least number of output data arrays
!>          4) The selected node is removed from the tree, the tree depths of the remaining data arrays are updated.
!>          5) if there is a node left, go to step 1.


!> \authors Stephan Thober
!> \date Oct 2019
#include "flogging.h"
module mo_mpr_reorder_data_array

  use mo_kind, only : i4
  use mo_mpr_constants, only : maxNoDataArrays, maxNameLength, noname
  use mo_mpr_data_array, only : DataArray
  use mo_constants, only: nodata_i4
  use flogging

  
  implicit none

  private

  public :: reorder_data_arrays

  ! --------------------------------------------------------------------------------------

contains

  subroutine reorder_data_arrays(MPR_DATA_ARRAYS)

    type(DataArray), intent(inout) :: MPR_DATA_ARRAYS(:) ! (maxNoDataArrays)
    
    ! variables for tree structure
    logical, allocatable :: mask(:)
    integer(i4) :: i, n
    integer(i4) :: sel_ind(1) ! selected index
    integer(i4) :: ind(1) ! helper index for copying
    integer(i4) :: n_tree ! number of elements in tree
    integer(i4) , allocatable :: tree_depth(:)
    integer(i4) , allocatable :: tree_order(:)
    character(maxNameLength) :: tree_inputFieldNames(maxNoDataArrays, 100)
    character(maxNameLength) :: tree_outputFieldNames(maxNoDataArrays, 100)
    character(maxNameLength) :: tree_deleteFieldNames(maxNoDataArrays, 100)
    character(maxNameLength), allocatable :: tree_name(:), tn(:), currentDeleteFieldNames(:)
    type(DataArray), allocatable :: DataArray_tmp(:)

    allocate(tree_depth(maxNoDataArrays))
    allocate(tree_order(maxNoDataArrays))
    allocate(tree_name(maxNoDataArrays))
    allocate(tn(maxNoDataArrays))
    
    ! initialize variables
    tree_name = noname
    tree_inputFieldNames = noname
    tree_outputFieldNames = noname
    tree_deleteFieldNames = noname
    tree_depth = nodata_i4
    tree_order = nodata_i4
    
    ! resort data arrays to minimize memory requirements
    do i = 1, maxNoDataArrays
      ! stop if no name is initialized
      if (.not. MPR_DATA_ARRAYS(i)%is_initialized) cycle
    
      tree_name(i) = MPR_DATA_ARRAYS(i)%name
      ! calculate initial tree depth of data array
      if (.not. MPR_DATA_ARRAYS(i)%transferHelper%doInitfromFile) then
        tree_inputFieldNames(i, :size(MPR_DATA_ARRAYS(i)%transferHelper%inputFieldNames)) = MPR_DATA_ARRAYS(i)%transferHelper%inputFieldNames
      end if
      call get_tree_depth(i, &
          MPR_DATA_ARRAYS(:)%name, &
          tree_inputFieldNames(i,:), &
          tree_outputFieldNames, &
          tree_depth)
    end do
    ! finished initialization

    ! remove nodata from tree
    mask = (tree_name /= noname)
    n_tree = count(mask)
    tree_name = pack(tree_name, mask)
    tree_depth = pack(tree_depth, mask)
    tree_order = pack(tree_order, mask)

    ! copy tree data that can be changed to stor
    tn = tree_name

    ! select next index
    do i = 1, n_tree
      log_trace(*) 'reorder_data_arrays: <<<<<<<<<<<<<<<<<'
      log_trace(*) 'reorder_data_arrays: iter: ', i

      ! print input tree
      do n = 1, size(tn)
        ! stop if no name is initialized
        if (tn(n) .eq. noname) cycle
        log_subtrace(*) 'reorder_data_arrays: scanned name: ', trim(tn(n)), ' with depth ', tree_depth(n), &
                ' and valid output fields: ', count(tree_outputFieldNames(n, :) .ne. noname)
      end do

      sel_ind = select_node(tn, tree_depth, tree_inputFieldNames, tree_outputFieldNames)
      log_trace(*) "reorder_data_arrays: selected id:   ", sel_ind(1)
      log_trace(*) "reorder_data_arrays: selected name: ", trim(tree_name(sel_ind(1))), tree_depth(sel_ind(1))
      
      ! set order
      tree_order(sel_ind(1)) = i

      ! remove item from tree by setting name to none and update tree
      tn(sel_ind(1)) = noname
      tree_outputFieldNames = noname
      do n = 1, n_tree
        call get_tree_depth(n, &
            tn, &
            tree_inputFieldNames(n, :), &
            tree_outputFieldNames, &
            tree_depth)
      end do

      ! determine items to delete
      call delete_arrays(sel_ind(1), n_tree, tree_name, tn, tree_inputFieldNames, tree_deleteFieldNames)
    end do

    if (n_tree .ne. count(tree_deleteFieldNames .ne. noname)) then
      log_error(*) 'size of tree: ', n_tree
      log_error(*) 'number of deleted items: ', count(tree_deleteFieldNames .ne. noname)
      log_error(*) 'ERROR*** number of deleted items is not equal to the number of items in data array'
      stop 1
    end if

    ! actually reorder DataArray and add delete FieldNames
    DataArray_tmp = MPR_DATA_ARRAYS
    do i = 1, n_tree
      call MPR_DATA_ARRAYS(i)%reset()
      ind = minloc(abs(tree_order - i))
      MPR_DATA_ARRAYS(i) = DataArray_tmp(ind(1))
      currentDeleteFieldNames = pack(tree_deleteFieldNames(ind(1), :), tree_deleteFieldNames(ind(1), :) /= noname)
      call MPR_DATA_ARRAYS(i)%add_delete_items(currentDeleteFieldNames)
      ! end info print out
      log_debug(*) 'reorder_data_arrays: moved DataArray "', trim(tree_name(ind(1))), '" from index ', ind(1), &
              ' to new index ', i, '. After handling of DataArray, ', size(currentDeleteFieldNames), &
              ' DataArrays are removed (parents and/or itself).'
      deallocate(currentDeleteFieldNames)
    end do

    ! free memory
    deallocate(mask)
    deallocate(tree_depth)
    deallocate(tree_order)
    deallocate(tree_name)
    deallocate(tn)
    deallocate(DataArray_tmp)
    
  end subroutine reorder_data_arrays

  subroutine get_tree_depth(ind, arrayNames, inputFieldNames, outputFieldNames, tree_depth)
  !< calculates the tree depth for each node and creates the outputFieldNames,
  !< i.e. those data arrays where the current data array is used

    integer(i4), intent(in) :: ind
    character(maxNameLength), intent(in) :: arrayNames(:)
    character(maxNameLength), intent(in) :: inputFieldNames(:)
    character(maxNameLength), intent(inout) :: outputFieldNames(:, :)
    integer(i4), intent(out) :: tree_depth(:)

    logical, allocatable :: input_mask(:)
    integer(i4) :: tree_depth_tmp(1)
    integer(i4) :: i, j, n, indOut(1)
    integer(i4) :: int_indexes(size(arrayNames, dim=1))

    ! get maximum tree depth of parents
    tree_depth(ind) = 0_i4
    do i = 1, size(int_indexes, dim=1)
      int_indexes(i) = i
    end do

    n = count(inputFieldNames /= noname)
    allocate(input_mask(size(arrayNames)))
    
    do i = 1, n
      log_subtrace(*) 'inputFieldNames: ', trim(inputFieldNames(i)), i, n
      !
      ! check whether inputfieldnames are within the tree, otherwise they were processed already
      input_mask = (arrayNames .eq. inputFieldNames(i))
      if (any(input_mask)) then
        ! add Name to OutputFieldName
        indOut = pack(int_indexes, input_mask)
        do j = 1, size(OutputFieldNames, dim=2)
          if (outputFieldNames(indOut(1), j) .eq. noname) then
            outputFieldNames(indOut(1), j) = arrayNames(ind)
            exit
          end if
        end do
        !
        ! store maximum tree depth of parents
        tree_depth_tmp = pack(tree_depth, input_mask)
        if (tree_depth_tmp(1) .gt. tree_depth(ind)) tree_depth(ind) = tree_depth_tmp(1)
      end if
    end do
    ! tree depth must be by one number higher than maximum of parents
    if (arrayNames(ind) /= noname) tree_depth(ind) = tree_depth(ind) + 1_i4
    
    deallocate(input_mask)
    
  end subroutine get_tree_depth

  function select_node(tn, td, tree_inputFieldNames, tree_outputFieldNames)
  !< The algorithm to select a node works like this:
  !< 1) algorithm always picks element with depth one.
  !< 2) if there are multiple, it selects the one where most previously calculated data arrays can be deleted
  !< 3) if nothing can be deleted, it picks the one with the least number of output data arrays
  !< 4) The selected node is removed from the tree, the tree depths of the remaining data arrays are updated.
  !< 5) if there is a node left, go to step 1.

    character(maxNameLength), intent(in) :: tn(:)
    integer(i4), intent(in) :: td(:)
    character(maxNameLength), intent(in) :: tree_inputFieldNames(maxNoDataArrays, 100)
    character(maxNameLength), intent(in) :: tree_outputFieldNames(maxNoDataArrays, 100)
    integer(i4) :: select_node(1)

    logical :: delete_item
    logical, allocatable :: mask_md(:)
    integer(i4) :: td_local(size(td, dim=1))
    integer(i4) :: j, k, l
    integer(i4) :: n_tree ! number of elements in tree
    integer(i4) , allocatable :: count_deletes(:), count_outputs(:)

    ! always select node with smallest depth
    n_tree = size(tn, dim=1)
    td_local = merge(huge(i4), td, td .eq. 0_i4)
    if (minval(td_local) .ne. 1) then
      log_error(*) "ERROR, minimum value of tree depth is not one, but: ", minval(td_local)
      stop 1
    end if
    mask_md = (td_local .eq. 1_i4)

    if (count(mask_md) .gt. 1) then

      allocate(count_outputs(size(mask_md)))
      allocate(count_deletes(size(mask_md)))
      count_outputs = 0_i4
      count_deletes = 0_i4

      do j = 1, n_tree
        if (.not. mask_md(j)) then
          count_outputs(j) = huge(i4)
          cycle
        end if

        ! count number of output fields, the lower -> the better
        count_outputs(j) = count(tree_outputFieldNames(j, :) .ne. noname)
      end do

      do j = 1, n_tree
        if (.not. mask_md(j)) cycle

        ! count if other arrays could be deleted if this array is selected
        if (count_outputs(j) .eq. 0) count_deletes(j) = count_deletes(j) + 1_i4
        do k = 1, n_tree
          if (tree_inputFieldNames(j, k) .eq. noname) exit
          
          delete_item = .true.
          do l = 1, n_tree
            if (tn(l) .eq. noname) cycle
            if (l .eq. j) cycle
            if (any(tree_inputFieldNames(l, :) .eq. tree_inputFieldNames(j, k))) delete_item = .False.
          end do
          if (delete_item) count_deletes(j) = count_deletes(j) + 1_i4
        end do
      end do
      log_trace(*) 'count_outputs: ', pack(count_outputs, mask_md)
      log_trace(*) 'count_deletes: ', count_deletes
      
      ! select maximum number of deletes
      if (maxval(count_deletes) .gt. 0_i4) then
        ! select maximum number of deletes
        log_trace(*) "selected maximum number of deletes: ", count(count_deletes .eq. maxval(count_deletes))
        select_node = maxloc(count_deletes)
      else
        log_trace(*) "selected minimum number of outputs"
        select_node = minloc(count_outputs)
      end if

      deallocate(count_outputs, count_deletes)
    else
      ! set an id
      log_trace(*) "selected minimum tree depth"
      select_node = minloc(td_local)
    end if

    ! check that tree depth is 1
    ! otherwise, it is not guaranteed that all input names are available
    if (.not. td(select_node(1)) .eq. 1_i4) then
      log_error(*) 'ERROR*** tree depth is not 1, instead: ', td(select_node(1))
      stop 1
    end if

    ! free memory
    deallocate(mask_md)
    
  end function select_node

  subroutine delete_arrays(sel_ind, n_tree, tree_name, tn, tree_inputFieldNames, tree_deleteFieldNames)
  !< checks whether any of the already calculated data arrays are still needed
  !< if not, they can be deleted

    integer(i4), intent(in) :: sel_ind
    integer(i4), intent(in) :: n_tree
    character(maxNameLength), intent(in) :: tree_name(:)
    character(maxNameLength), intent(in) :: tn(:)
    character(maxNameLength), intent(in) :: tree_inputFieldNames(maxNoDataArrays, 100)
    character(maxNameLength), intent(out) :: tree_deleteFieldNames(maxNoDataArrays, 100)

    integer(i4) :: j, k
    logical :: delete_item
    
    do j = 1, n_tree
      ! j-th data array has not been calculated -> cannot be deleted
      if (tn(j) .ne. noname) cycle
      ! j-th data array has been already deleted -> cannot be deleted
      if (any(tree_deleteFieldNames .eq. tree_name(j))) cycle
      !
      ! check whether j-th data array should be deleted
      delete_item = .true.
      ! check whether it has been calculated and is needed in the remaining nodes
      do k = 1, n_tree
        if (tn(k) .eq. noname) cycle
        if (any(tree_inputFieldNames(k, :) .eq. tree_name(j))) delete_item = .false.
      end do
      ! if j-th item should be deleted, append it to deleteFieldNames
      if (delete_item) then
        do k = 1, size(tree_deleteFieldNames, 2)
          if (tree_deleteFieldNames(sel_ind, k) .eq. noname) then
            tree_deleteFieldNames(sel_ind, k) = tree_name(j)
            exit
          end if
        end do
      end if
    end do
  end subroutine delete_arrays

end module mo_mpr_reorder_data_array
