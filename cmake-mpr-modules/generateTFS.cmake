# generate transfer function fortran code with python on the fly
# file paths will be converted to absolute paths with this:
set(GENERATE_TFS_CONFIG "" CACHE FILEPATH "path to config file for MPR")
set(GENERATE_TFS_PARAM "" CACHE FILEPATH "path to config file with extra parameters for MPR")
if(NOT (GENERATE_TFS_CONFIG STREQUAL ""))
  message(STATUS "MPR: generate TFs from config '${GENERATE_TFS_CONFIG}'")
  # check for python
  find_package(Python3 COMPONENTS Interpreter REQUIRED)
  if (DEFINED ENV{PYTHON})
    set(Python3_EXECUTABLE $ENV{PYTHON})
    message(STATUS "Overwriting Python3_EXECUTABLE to " ${Python3_EXECUTABLE} " as set in env var 'PYTHON'")
  endif()
  # check for f90nml
  execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import f90nml"
    RESULT_VARIABLE EXIT_CODE
    OUTPUT_QUIET
  )
  if (NOT ${EXIT_CODE} EQUAL 0)
    message(
      FATAL_ERROR
      "The \"f90nml\" Python3 package is not installed. Please install it using the following command: \"pip3 install f90nml\"."
    )
  endif()
  # check for additional param file
  if(NOT (GENERATE_TFS_PARAM STREQUAL ""))
    message(STATUS "MPR: generate TFs with additional parameter from '${GENERATE_TFS_PARAM}'")
    # run it
    execute_process(
      COMMAND ${Python3_EXECUTABLE} -m src_python.pre_proc.update_tfs_in_fortran_source -c ${GENERATE_TFS_CONFIG} -p ${GENERATE_TFS_PARAM}
      # currently in src/
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/..)
  else()
    message(STATUS "MPR: generate TFs with no additional parameter")
    # run it
    execute_process(
      COMMAND ${Python3_EXECUTABLE} -m src_python.pre_proc.update_tfs_in_fortran_source -c ${GENERATE_TFS_CONFIG}
      # currently in src/
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/..)
  endif()
else()
  message(STATUS "MPR: no generation of TFs from config")
endif()
