# !/bin/bash

cd ../cfilt
clang-format -i *.c *.h

cd ../tests
clang-format -i *.c *.h

cd ../scripts
