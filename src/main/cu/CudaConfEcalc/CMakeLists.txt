cmake_minimum_required(VERSION 3.10)

project(CudaConfEcalc VERSION 0.1 DESCRIPTION "CUDA-based conformation energy calculator for OSPREY")

# NOTE: if CMake complains about not being able to find nvcc, try something like this:
#set(CMAKE_CUDA_COMPILER /opt/cuda-10.2/bin/nvcc)

# turn on CUDA support
enable_language(CUDA)

# *** Don't Try To Compile This Code With Cuda Toolkit v11 or Later Unless You Know What You're Doing !! ***
# There must be some kind of regression in nvcc v11+, because the compiled code can no longer
# meet the launch bounds on kernels, whereas the v10.2 compiler has no problems.
# It looks like v11+ uses too many registers and goes waaaaay over the register limits imposed by the launch bounds.
# The only known working compiler version is 10.2.
# Use other versions at your own risk.

set(CMAKE_CUDA_STANDARD 14)

# let CMake create defines we can use in the code
configure_file(config.h.in config.h)

# define the source files to compile
set(CXXFiles
        confecalc.cu
        cuda.cu
        motions/dihedral.cu
        motions/transrot.cu
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
        CUDA_ARCHITECTURES false # newer versions of CMake require this, even though we don't use it
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
        CUDA_ARCHITECTURES false # newer versions of CMake require this, even though we don't use it
        )
