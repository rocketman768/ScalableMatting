#!/bin/bash

for im in 01 02 03 04 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27
do
   convert "GT${im}.png" "${im}_gt.pgm"
done

