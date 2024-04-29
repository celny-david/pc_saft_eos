# DYNAMIC/SHARED library with the test of execution for PC_SAFT
# resource handling is build into the program to run from the make directory
# or primarily run resources near the library

project(LIB_DYNAMIC)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}) # this work for shared library
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})

add_library(pcsaft_dynamic_lib SHARED "${MOD_INTERFACE}" "${MOD_SOURCES}" "${SUB_SOURCES}")
set_target_properties(pcsaft_dynamic_lib PROPERTIES OUTPUT_NAME "pcsaft") # the lib prefix is added automatically
#message("lib_out  : " "${PROJECT_SOURCE_DIR}" )

add_executable(pcsaft_dynamic_exec ${CMAKE_SOURCE_DIR}/test/test_fortran_library_dynamic.f90)
target_link_libraries(pcsaft_dynamic_exec pcsaft_dynamic_lib )
