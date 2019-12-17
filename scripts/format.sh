# !/bin/bash

CUR=$(pwd)
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

cd "$DIR"

cd ../cfilt
pwd
clang-format -i *.c *.h

cd ../tests
pwd
clang-format -i *.c *.h
