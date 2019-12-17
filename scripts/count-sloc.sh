# !/bin/bash

CUR=$(pwd)
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

cd "$DIR"
cd ..

rm -rf bin build
find -name "*.sh" -or -name "*.c" -or -name "*.h" -or -name "*.py" | xargs wc -l
