# !/bin/bash

cd ..
rm -rf bin build
find -name "*.sh" -or -name "*.c" -or -name "*.h" | xargs wc -l
