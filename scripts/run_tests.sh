# !/bin/bash

CUR=$(pwd)
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

cd "$DIR"
cd ..
bin/test_gh       > tests/test_gh.csv
bin/test_kalman1d > tests/test_kalman1d.csv
bin/test_kalman   > tests/test_kalman.csv

cd "$CUR"
