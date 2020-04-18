cmake_minimum_required(VERSION 3.10)

project(ConfEcalc VERSION 0.1 DESCRIPTION "native conformation energy calculator for OSPREY")

set(CMAKE_C_STANDARD 99)

# let CMake create defines we can use in the code
configure_file(config.h.in config.h)

# define the source files to compile
set(CFiles
        confEcalc.c
        confSpace.c
        math.c
        )

# get paths into the parent Osprey project
set(OspreyPath ../../../..)
set(JnaPlatform linux-x86-64) # TODO: get the JNA names for each platform
set(OspreyBuildPath ${OspreyPath}/build/resources/main/${JnaPlatform})

# tell CMake to write the libs into Osprey's resources dir
set(LIBRARY_OUTPUT_PATH ${OspreyPath}/resources/${JnaPlatform})

# make targets for each value of REAL
foreach(Real f64 f32)

    # make a library target
    add_library(ConfEcalc_${Real} SHARED ${CFiles})
    target_compile_definitions(ConfEcalc_${Real} PRIVATE REAL=${Real})
    target_include_directories(ConfEcalc_${Real} PUBLIC "${PROJECT_BINARY_DIR}")

    # add another target to copy the libraries from the resources path to the build path
    # so we can run java with newly-compiled libraries without having to run gradle to process the resources
    add_custom_target(ConfEcalc_${Real}_CopyLibs)
    add_dependencies(ConfEcalc_${Real}_CopyLibs ConfEcalc_${Real})
    add_custom_command(
            TARGET ConfEcalc_${Real}_CopyLibs POST_BUILD
            COMMAND cp ${LIBRARY_OUTPUT_PATH}/$<TARGET_FILE_NAME:ConfEcalc_${Real}> ${OspreyBuildPath}
    )

endforeach(Real)
