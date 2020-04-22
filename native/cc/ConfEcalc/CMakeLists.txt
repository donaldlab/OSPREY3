cmake_minimum_required(VERSION 3.10)

project(ConfEcalc VERSION 0.1 DESCRIPTION "native conformation energy calculator for OSPREY")

set(CMAKE_CXX_STANDARD 17)

# let CMake create defines we can use in the code
configure_file(config.h.in config.h)

# define the source files to compile
set(CXXFiles
        confecalc.cc
        )

SET(CMAKE_CXX_FLAGS -pedantic-errors)

# get paths into the parent Osprey project
set(OspreyPath ../../../..)
set(JnaPlatform linux-x86-64) # TODO: get the JNA names for each platform
set(OspreyBuildPath ${OspreyPath}/build/resources/main/${JnaPlatform})

# tell CMake to write the libs into Osprey's resources dir
set(LIBRARY_OUTPUT_PATH ${OspreyPath}/resources/${JnaPlatform})

# make a library target
add_library(ConfEcalc SHARED ${CXXFiles})
target_include_directories(ConfEcalc PUBLIC "${PROJECT_BINARY_DIR}")

# add another target to copy the libraries from the resources path to the build path
# so we can run java with newly-compiled libraries without having to run gradle to process the resources
add_custom_target(ConfEcalc_CopyLibs)
add_dependencies(ConfEcalc_CopyLibs ConfEcalc)
add_custom_command(
        TARGET ConfEcalc_CopyLibs POST_BUILD
        COMMAND cp ${LIBRARY_OUTPUT_PATH}/$<TARGET_FILE_NAME:ConfEcalc> ${OspreyBuildPath}
)

