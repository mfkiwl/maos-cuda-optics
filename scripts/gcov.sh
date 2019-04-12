#!/bin/sh
#To assess coverage, do the following
#1. compile maos with --enable-gcov
#2. Run maos
#3. run this script at the root of the building directory
#4. view out/index.html
lcov --capture --directory . - coverage.info || exit
lcov --remove coverage.info '/usr/*' '/opt/*' '/tmp/*' -o coverage2.info || exit
genhtml coverage2.info - out
