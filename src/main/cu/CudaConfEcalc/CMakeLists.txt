cmake_minimum_required(VERSION 3.10)

project(CudaConfEcalc VERSION 0.1 DESCRIPTION "CUDA-based conformation energy calculator for OSPREY")

# turn on CUDA support
enable_language(CUDA)

set(CMAKE_CUDA_STANDARD 14)

# let CMake create defines we can use in the code
configure_file(config.h.in config.h)

# define the source files to compile
set(CXXFiles
        confecalc.cu
        motions.cu
        cuda.cu
        )

SET(CMAKE_CXX_FLAGS -pedantic-errors)

# get paths into the parent Osprey project
set(OspreyPath ${CMAKE_SOURCE_DIR}/../../../..)
set(JnaPlatform linux-x86-64) # TODO: get the JNA names for each platform
set(OspreyBuildPath ${OspreyPath}/build/resources/main/${JnaPlatform})

# tell CMake to write the libs into Osprey's resources dir
set(LIBRARY_OUTPUT_PATH ${OspreyPath}/src/main/resources/${JnaPlatform})

# set the CUDA options to build for all modern GPU architectures
# see compatibilty guide:
# https://docs.nvidia.com/cuda/turing-compatibility-guide/index.html
set(CMAKE_CUDA_FLAGS "\
        -Xcudafe --display_error_number \
        --resource-usage \
        --generate-line-info \
        -gencode=arch=compute_50,code=sm_50 \
        -gencode=arch=compute_52,code=sm_52 \
        -gencode=arch=compute_60,code=sm_60 \
        -gencode=arch=compute_61,code=sm_61 \
        -gencode=arch=compute_70,code=sm_70 \
        -gencode=arch=compute_75,code=sm_75 \
        -gencode=arch=compute_75,code=compute_75 \
        ")

# make a library target, for calling from osprey
add_library(CudaConfEcalc SHARED ${CXXFiles} api.cu)
target_include_directories(CudaConfEcalc PUBLIC "${PROJECT_BINARY_DIR}")
set_target_properties(CudaConfEcalc PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        POSITION_INDEPENDENT_CODE ON
        )

# add another target to copy the libraries from the resources path to the build path
# so we can run java with newly-compiled libraries without having to run gradle to process the resources
add_custom_target(CudaConfEcalc_CopyLibs)
add_dependencies(CudaConfEcalc_CopyLibs CudaConfEcalc)
add_custom_command(
        TARGET CudaConfEcalc_CopyLibs POST_BUILD
        COMMAND cp ${LIBRARY_OUTPUT_PATH}/$<TARGET_FILE_NAME:CudaConfEcalc> ${OspreyBuildPath}
)

# make an executable target, for running job dumps without osprey
add_executable(rundump ${CXXFiles} rundump.cu)
target_include_directories(rundump PUBLIC "${PROJECT_BINARY_DIR}")
set_target_properties(rundump PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        POSITION_INDEPENDENT_CODE ON
        )
