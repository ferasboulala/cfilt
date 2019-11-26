# !/bin/bash

cd ../cfilt
pwd
clang-format -i *.c *.h

cd ../tests
pwd
clang-format -i *.c

cd ../scripts
pwd
