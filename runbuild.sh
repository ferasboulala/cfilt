#!/bin/bash

mkdir -p build bin
cd build
cmake ..
make -j4
