# !/bin/bash

CUR=$(pwd)
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

cd "$DIR"

cd ../cfilt
pwd
clang-format -i *.c *.h

cd ../tests
pwd
find -name "*.h" -or -name "*.c" | xargs clang-format -i

cd ../examples
pwd
find -name "*.h" -or -name "*.c" | xargs clang-format -i
