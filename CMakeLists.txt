cmake_minimum_required(VERSION 3.14)
# check version file
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake-modules)
# get version and date from files (version.txt and version_date.txt)
include(version)
get_version(MPR_VER MPR_VER_DEV MPR_DATE)

# create the project
project(MPR
        VERSION ${MPR_VER}
        DESCRIPTION "Multiscale Parameter Regionalization"
        HOMEPAGE_URL "https://www.ufz.de/index.php?en=40126"
        LANGUAGES Fortran)

# add full version and date to pre-processor flags (qoutes need in before hand)
add_compile_definitions(MPRVERSION='${MPR_VER_DEV}' MPRDATE='${MPR_DATE}')

# common (default) options to cmake
# -DCMAKE_BUILD_TYPE=Release     - compile in debug or release mode
# -DBUILD_TESTING=ON             - whether to compile tests
# -DCMAKE_VERBOSE_MAKEFILE=OFF   - see all the commands
# -DCMAKE_BUILD_MODULE_SYSTEM_INDEPENDENT=OFF   - build the library INDEPENDENT of the module system,
# -DCMAKE_FIND_LIBRARY_CUSTOM_LIB_SUFFIX=64     - find extra library paths on some systems (eve.nag62)
# -DGENERATE_TFS_CONFIG="mpr_mhm.nml"  - "path to config file for MPR"
# -DGENERATE_TFS_PARAM="mpr_global_parameter_mhm.nml" - "path to config file with extra parameters for MPR"

# automatically enables testing
include(CTest)

# The variable "CMAKE_BUILD_MODULE_SYSTEM_INDEPENDENT" can be set before executing cmake via a cache command:
# $cmake -DCMAKE_BUILD_MODULE_SYSTEM_INDEPENDENT:STRING=ON ..
# or cache file:
# $cmake -C ../CMakeCacheFiles/eve ..
# or after executing CMake editing the CMakeCache.txt, preferably with a corresponding cmake editor i.e ccmake
option(CMAKE_BUILD_MODULE_SYSTEM_INDEPENDENT "build the module INDEPENDENT of the module system, so the build in the build tree works even after a module purge")

# same with coverage
option(CMAKE_WITH_COVERAGE "build the module with gcov coverage support")

# library module specific settings
add_subdirectory(./src)
add_subdirectory(./src_python)

# this command is able to create dependencies, compile and add the sources in the right order
set(EXECUTABLE_NAME "${CMAKE_PROJECT_NAME}_runner")
add_executable(${EXECUTABLE_NAME} ./src/main.F90)
target_link_libraries(${EXECUTABLE_NAME} PRIVATE mpr)

set_target_properties ( ${EXECUTABLE_NAME}
  PROPERTIES
    OUTPUT_NAME "${CMAKE_PROJECT_NAME}"
    VERSION "${MPR_VERSION}"
)
if (CMAKE_BUILD_MODULE_SYSTEM_INDEPENDENT)
  set_target_properties ( ${EXECUTABLE_NAME}
    PROPERTIES
      INSTALL_RPATH "${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}"
      BUILD_WITH_INSTALL_RPATH ON
  )
endif()
