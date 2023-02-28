#!/bin/zsh


mkdir $1
for x in 1 2 4 6 8 10 13 16 19 20 30 40 50 60 70 80 90 100; do touch "./$1/t1_$2_$x.csv"; done

touch "./$1/t2_$2_p.csv"