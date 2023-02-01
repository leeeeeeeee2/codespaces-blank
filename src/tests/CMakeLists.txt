enable_testing()

set (LIB_NAME mpr)
find_package(PFUNIT QUIET)

if(PFUNIT_FOUND)
	# find if there is an OpenMP setup on the system and if so, set corresponding variables
	find_package(OpenMP)
	if (NOT ${OpenMP_Fortran_FOUND})
		message(FATAL_ERROR "OpenMP required but not found")
	endif()

	file ( COPY ./files
	DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

  set (pfunit_suffix ".pf")
  file( GLOB testsourcefiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *${pfunit_suffix} )
  # set (testnames "")
  foreach( testsourcefile ${testsourcefiles} )
    # I used a simple string replace, to cut off .pf
    string( REPLACE ${pfunit_suffix} "" testname ${testsourcefile} )
    add_pfunit_ctest (${testname}
      TEST_SOURCES ${testsourcefile}
      LINK_LIBRARIES ${LIB_NAME}
      )
    list(APPEND testnames ${testname})
  endforeach()

  if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    if (CMAKE_WITH_COVERAGE)
      include(CodeCoverage)
      APPEND_COVERAGE_COMPILER_FLAGS()
      # add_custom_target(test_runner COMMAND ${testnames})
      SETUP_TARGET_FOR_COVERAGE_LCOV(NAME MPR_coverage_CI
        EXECUTABLE ${testnames}
        # BASE_DIRECTORY "${PROJECT_SOURCE_DIR}/src"
        # DEPENDENCIES {LIB_NAME}
        GENHTML_ARGS -t "${LIB_NAME} coverage" --html-prolog ${CMAKE_CURRENT_SOURCE_DIR}/../../doc/html_files/cov_header.prolog)
    endif()
  endif()

else()
    message(STATUS "No test program found, install PFUNIT")
endif()
