#!/bin/bash

for i in PM-*/*/
do
cp ktune.sh $i
done

for t in PM-*/*
do
cd $t
sbatch *sh
cd ../..
done
