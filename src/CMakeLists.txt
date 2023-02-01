include(FetchContent)

# generate transfer function fortran code with python on the fly
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../cmake-mpr-modules)
include(generateTFS)

set (LIB_NAME mpr)
file (GLOB sources ./mo_mpr_*.f90 ./mo_mpr_*.F90)

option(BUILD_MPR_LIB_SHARED "Build MPR library as shared." OFF)
if(BUILD_MPR_LIB_SHARED)
  add_library(${LIB_NAME} SHARED ${sources})
else()
  add_library(${LIB_NAME} STATIC ${sources})
endif()
target_include_directories(${LIB_NAME} PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
# target_include_directories(${LIB_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

include(../cmake/cmake-modules/CPM.cmake)

set (TARGET_NAME forces)
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../${TARGET_NAME}")
  message(STATUS "${LIB_NAME}: found local ${TARGET_NAME} directory")
  set(CPM_forces_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/../${TARGET_NAME}" CACHE PATH "Local source path for ${TARGET_NAME}.")
else()
  set(CPM_forces_SOURCE "" CACHE PATH "Local source path for ${TARGET_NAME}.")
endif()
CPMAddPackage("https://git.ufz.de/chs/forces.git#main")
if(BUILD_MPR_LIB_SHARED)
  set_property(TARGET ${TARGET_NAME} PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()
target_link_libraries(${LIB_NAME} PUBLIC ${TARGET_NAME})

set (TARGET_NAME flogging)
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../${TARGET_NAME}")
  message(STATUS "${LIB_NAME}: found local ${TARGET_NAME} directory")
  set(CPM_forces_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/../${TARGET_NAME}" CACHE PATH "Local source path for ${TARGET_NAME}.")
else()
  set(CPM_forces_SOURCE "" CACHE PATH "Local source path for ${TARGET_NAME}.")
endif()
CPMAddPackage("https://github.com/MuellerSeb/flogging.git#pre_proc_fix")
set (${FLOGGING_NAME}_BUILD_TESTING OFF)
if(BUILD_MPR_LIB_SHARED)
  set_property(TARGET ${TARGET_NAME} PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()
target_link_libraries(${LIB_NAME} PUBLIC ${TARGET_NAME})


# by setting compile options and definitions PUBLIC, they are also used by
# programms linking agains it (mhm exe in this case)
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  target_compile_definitions(${LIB_NAME} PRIVATE GFORTRAN)
  target_compile_options(${LIB_NAME} PUBLIC
    -ffixed-line-length-none
    $<$<CONFIG:DEBUG>:-Og -g -Wall -Wextra -fimplicit-none -fbacktrace -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=snan>
    $<$<CONFIG:RELEASE>:-Ofast>
    $<$<BOOL:${CMAKE_WITH_COVERAGE}>:-g --coverage>
  )
  if(CMAKE_WITH_COVERAGE)
    target_link_libraries(${LIB_NAME} PUBLIC gcov)
  endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  # https://discourse.cmake.org/t/preserving-options-with-spaces-in-add-compile-options/1551/2
  target_compile_definitions(${LIB_NAME} PRIVATE INTEL)
  target_compile_options(${LIB_NAME} PUBLIC
    -nofixed "SHELL:-assume byterecl" "SHELL:-fp-model source" -m64 "SHELL:-assume realloc-lhs"
    $<$<CONFIG:DEBUG>:-g "SHELL:-warn all" "SHELL:-check all" -debug extended -traceback -fp-stack-check -O0 -fstack-protector-all -fstack-security-check -fpe0>
    $<$<CONFIG:RELEASE>:-O3 -qoverride-limits>
  )
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "NAG")
  target_compile_definitions(${LIB_NAME} PRIVATE NAG)
  target_compile_options(${LIB_NAME} PUBLIC
    -colour -unsharedf95 -ideclient -fpp
    $<$<CONFIG:DEBUG>:-g -nan -O0 -C=all -strict95 -ieee=stop>
    $<$<CONFIG:RELEASE>:-O4 -ieee=full>
  )
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
  target_compile_definitions(${LIB_NAME} PRIVATE NAG)
  target_compile_options(${LIB_NAME} PUBLIC
    -Mfree -mcmodel=medium
    $<$<CONFIG:DEBUG>:-C -c -g -traceback -O0>
    $<$<CONFIG:RELEASE>:-fast>
  )
endif()

# ifort and gfortran need the flag -cpp to interpret definitions like -DMRM2MHM
# the nag compiler is not able to interpret the flag -cpp but can interpret these definitions anyway
# so we check whether the compiler is able to use the flag -cpp
# for that we need the module CheckFortranCompilerFlag
include(CheckFortranCompilerFlag)
CHECK_Fortran_COMPILER_FLAG("-cpp" CPP_FLAG)
# if the flag exists, we add it to the compilation flags
if (CPP_FLAG)
	target_compile_options(${LIB_NAME} PUBLIC "-cpp")
endif()

if((CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME OR ${LIB_NAME}_BUILD_TESTING) AND BUILD_TESTING)
  add_subdirectory(./tests)
endif()
