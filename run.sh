#!/bin/bash

SHIFT=4096
declare -a MATRIX_SIZES=(3 4 5 6 7 8 9 10 20 30 40 50 60 70 80 90 
                         100 200 300 500 1000 1500)

for i in "${MATRIX_SIZES[@]}"
do
    build/bin/hankelhacker "$i" "$SHIFT" >>results.txt
    echo >>results.txt
done
