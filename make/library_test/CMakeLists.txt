# DYNAMIC/SHARED library with the test of execution for PC_SAFT
# resource handling is build into the program to run from the make directory
# or primarily run resources near the library

project(LIB_TEST)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}) # this work for shared library
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})

set (PFUNIT_INSTALL_DIR ${CMAKE_BINARY_DIR}/pfunit/install)

# This assumes that you've already cloned the pFUnit source code, built and installed it (make install)
# to a suitable location on your computer. Adding search paths using the PATHS argument is often necessary.
# See https://cmake.org/cmake/help/v3.0/command/find_package.html for more information.
# Note that the pFUnit package will be OS, compiler and compiler version dependent as it distributes Fortran mod-files.
find_package("pFUnit" QUIET HINTS ${PFUNIT_INSTALL_DIR})

# Alternatively pFUnit could be built during configuration of the CMake build folder by commenting out the line above
# and uncommenting the following lines:
if(PFUNIT_FOUND)
 message ("pFUnit found")
else()
 set (PFUNIT_SRC_DIR ${CMAKE_BINARY_DIR}/pfunit/src)
 set (PFUNIT_DOWNLOAD_DIR ${CMAKE_BINARY_DIR}/pfunit/download)
 set (PFUNIT_BUILD_DIR ${CMAKE_BINARY_DIR}/pfunit/build)

 set (pFUnit_DIR ${PFUNIT_INSTALL_DIR})
 message ("Building pFUnit")
 include (${CMAKE_CURRENT_LIST_DIR}/build_pfunit_from_source.cmake)
 message ("Looking for pFUnit here: ${PFUNIT_INSTALL_DIR}")
 find_package ("pFUnit"
     CONFIG
     PATHS ${PFUNIT_INSTALL_DIR})
endif()

enable_testing()

# Now that we have found the pFUnit package we can define a test executbale.
# Here, it's assumed that we'd like to test a library defined in the LIB_NAME variable,
# e.g. set (LIB_NAME some_lib)
add_library(pcsaft_test_lib "${MOD_SOURCES}" "${SUB_SOURCES}" "${MOD_INTERFACE}")

#set (TEST_SOURCES2 file(RELATIVE_PATH ${TEST_SOURCES}))
file(RELATIVE_PATH TEST_SOURCES ${PROJECT_SOURCE_DIR} ${TEST_SOURCES})

#message("PROJECT_BINARY_DIR  : " ${PROJECT_BINARY_DIR}) # DEBUG for correct inlcude directory path check

# Add pFUnit test if any TEST_SOURCES is defined
if (NOT "${TEST_SOURCES}" STREQUAL "")

    add_pfunit_test(pcsaft_test
                    ${TEST_SOURCES}
                    ""
                    "")
    target_include_directories(pcsaft_test PRIVATE ${PROJECT_BINARY_DIR}/iclude)
    target_link_libraries (pcsaft_test
                           pcsaft_test_lib)
endif (NOT "${TEST_SOURCES}" STREQUAL "")
 