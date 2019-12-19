#!/bin/bash

kernels=(residueForcefield residueCcd)

for kernel in ${kernels[*]}
do
    nvcc -fatbin \
        -arch=sm_70 \
        -gencode=arch=compute_50,code=sm_50 \
        -gencode=arch=compute_52,code=sm_52 \
        -gencode=arch=compute_60,code=sm_60 \
        -gencode=arch=compute_61,code=sm_61 \
        -gencode=arch=compute_70,code=sm_70 \
        -gencode=arch=compute_75,code=sm_75 \
        -gencode=arch=compute_75,code=compute_75 \
        $kernel.cu -o $kernel.bin
done
