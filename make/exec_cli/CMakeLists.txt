# CLI execution variant of PC_SAFT
# resource handling is build into the program to run from the make directory
# or primarily run resources near the executable
project(CLI)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/make/exec_cli)

add_executable(pcsaft "${MOD_SOURCES}" "${SUB_SOURCES}" "${PROG_SOURCES}" )
