cmake_minimum_required(VERSION 3.10)

# make sure we're using Intel's compiler for this library
if(CMAKE_CXX_COMPILER STREQUAL "icc")
    message(STATUS "using Intel compiler")
else()
    message(FATAL_ERROR "\"${CMAKE_CXX_COMPILER}\" configured as compiler instead of Intel. try -DCMAKE_CXX_COMPILER=icc")
endif()

project(IntelConfEcalc VERSION 0.1 DESCRIPTION "native conformation energy calculator for OSPREY, for Intel platforms")

set(CMAKE_CXX_STANDARD 17)

# let CMake create defines we can use in the code
configure_file(config.h.in config.h)

# define the source files to compile
set(CXXFiles
        confecalc.cc
        )

# set Intel compiler flags
set(CMAKE_CXX_FLAGS "\
        -pedantic-errors \
        -static-intel \
        -qopt-report=5 \
        -qopt-report-phase=vec \
        ")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_DEBUG "-g")

# get paths into the parent Osprey project
set(OspreyPath ../../../..)
set(JnaPlatform linux-x86-64) # TODO: get the JNA names for each platform
set(OspreyBuildPath ${OspreyPath}/build/resources/main/${JnaPlatform})

# tell CMake to write the libs into Osprey's resources dir
set(LIBRARY_OUTPUT_PATH ${OspreyPath}/resources/${JnaPlatform})

# make a library target
add_library(IntelConfEcalc SHARED ${CXXFiles})
target_include_directories(IntelConfEcalc PUBLIC "${PROJECT_BINARY_DIR}")

# add another target to copy the libraries from the resources path to the build path
# so we can run java with newly-compiled libraries without having to run gradle to process the resources
add_custom_target(IntelConfEcalc_CopyLibs)
add_dependencies(IntelConfEcalc_CopyLibs IntelConfEcalc)
add_custom_command(
        TARGET IntelConfEcalc_CopyLibs POST_BUILD
        COMMAND cp ${LIBRARY_OUTPUT_PATH}/$<TARGET_FILE_NAME:IntelConfEcalc> ${OspreyBuildPath}
)

