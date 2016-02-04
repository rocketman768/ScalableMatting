#!/bin/bash

for im in 01 02 03 04 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27
do
   convert "../8Mpx/${im}.ppm" -resize 2000000@ "${im}.ppm"
   convert "../8Mpx/${im}_scribs.pgm" -resize 2000000@ "${im}_scribs.pgm"
   convert "../8Mpx/${im}_gt.pgm" -resize 2000000@ "${im}_gt.pgm"
done

