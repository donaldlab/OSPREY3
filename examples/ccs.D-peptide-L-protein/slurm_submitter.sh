#!/bin/bash

for file in $1/*/$1-*sh;
do
sbatch $file
done
