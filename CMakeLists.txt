# cmake file
cmake_minimum_required(VERSION 3.14)

project(PC_SAFT_GLOBAL
        VERSION 0.2
        DESCRIPTION "Implementation of the PC-SAFT equation of state by D. Celný"
        LANGUAGES Fortran)
enable_language(Fortran)
#set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f90)

# === the compiler flags ===
if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-ffree-line-length-none -std=f2008 -ffree-form")
    set(bounds "-fbounds-check")
    set(debug "-g -fbacktrace")
endif ()
if (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(dialect "-stand f08 -free -implicitnone")
    set(bounds "-check bounds")
endif ()
if (CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    set(dialect "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03")
    set(bounds "-C")
endif ()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${dialect} ${bounds} ${debug} -gdwarf-3")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")

# === sources listing ===
#file(GLOB_RECURSE COMMON_SOURCES src/*mod.f90 src/*_sub.f90)
if (WIN32)
    set (MOD_INTERFACE  ${CMAKE_SOURCE_DIR}/src/bind_interface_win_mod.f90)
endif (WIN32)
if (UNIX)
    set (MOD_INTERFACE  ${CMAKE_SOURCE_DIR}/src/bind_interface_mod.f90)
endif (UNIX)
set (MOD_SOURCES    ${CMAKE_SOURCE_DIR}/src/contrib_mod.f90
                    ${CMAKE_SOURCE_DIR}/src/control_mod.f90
                    ${CMAKE_SOURCE_DIR}/src/ideal_part_mod.f90
                    ${CMAKE_SOURCE_DIR}/src/interface_mod.f90
                    ${CMAKE_SOURCE_DIR}/src/param_list_mod.f90)
set (SUB_SOURCES    ${CMAKE_SOURCE_DIR}/src/a_disp_sub.f90
                    ${CMAKE_SOURCE_DIR}/src/a_hc_sub.f90
                    ${CMAKE_SOURCE_DIR}/src/a_hs_sub.f90
                    ${CMAKE_SOURCE_DIR}/src/z_disp_sub.f90
                    ${CMAKE_SOURCE_DIR}/src/z_hc_sub.f90
                    ${CMAKE_SOURCE_DIR}/src/z_hs_sub.f90
                    ${CMAKE_SOURCE_DIR}/src/z_dd_sub.f90
                    ${CMAKE_SOURCE_DIR}/src/z_dq_sub.f90
                    ${CMAKE_SOURCE_DIR}/src/z_qq_sub.f90
                    ${CMAKE_SOURCE_DIR}/src/z_asoc_sub.f90
                    ${CMAKE_SOURCE_DIR}/src/rdf_sub.f90
                    ${CMAKE_SOURCE_DIR}/src/num_diff_calc_sub.f90
                    ${CMAKE_SOURCE_DIR}/src/press_calc_sub.f90
                    ${CMAKE_SOURCE_DIR}/src/chem_pot_calc_sub.f90
                    ${CMAKE_SOURCE_DIR}/src/helmholtz_calc_sub.f90
                    ${CMAKE_SOURCE_DIR}/src/enthalpy_entropy_calc_sub.f90)
set (PROG_SOURCES   ${CMAKE_SOURCE_DIR}/src/pcsaft_main_prog.f90)

#message("common sources" "${COMMON_SOURCES}" )
#message("app sources" "${APP_SOURCES}" )

# === project path listing ===

add_subdirectory(make/exec_cli)
add_subdirectory(make/library_dynamic)
