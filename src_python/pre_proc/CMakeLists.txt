enable_testing()
# https://stackoverflow.com/questions/52494425/running-a-pytest-test-from-cmake-where-the-test-and-sources-are-in-different-fol
# https://www.scivision.dev/cmake-python-pytest/

find_package(Python3 COMPONENTS Interpreter)

if(PYTHON3_FOUND)
  # TODO: this is another dirty hack, as CMake has problems detecting the right Python version in lieu of
  # conda or virtualenv
  if (DEFINED ENV{PYTHON})
    set(Python3_EXECUTABLE $ENV{PYTHON})
    message(STATUS "Overwriting Python3_EXECUTABLE to " ${Python3_EXECUTABLE} " as set in env var 'PYTHON'")
  endif()
	add_test (NAME pytests
	  COMMAND ${Python3_EXECUTABLE} -m pytest # ${CMAKE_CURRENT_SOURCE_DIR}
	  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
	)
else()
    message(STATUS "No python/pytest environment found, testing skipped")
endif()
