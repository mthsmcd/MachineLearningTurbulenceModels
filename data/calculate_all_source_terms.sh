#!/bin/bash

sd_dirs="2200 2400 2600 2900 3200 3500"
ph_dirs="0p5 0p8 1p0 1p2 1p5"
fields="nut t tStar Rperp"

for dir in ${sd_dirs};
do

  cp -rv ./square-duct/"$dir"/0/Udns ./square-duct/"$dir"/0/U
  cp -rv ./square-duct/"$dir"/0/Rdns ./square-duct/"$dir"/0/R
  mv -v ./square-duct/"$dir"/0/nut ./square-duct/"$dir"/0/nutRans

  calculateNut -case ./square-duct/"$dir"

  calculateRFV -case ./square-duct/"$dir"

  calculateRFVperp -case ./square-duct/"$dir"

  calculateRperp -case ./square-duct/"$dir"

  for f in ${fields};
  do

    mv -v ./square-duct/"$dir"/0/"$f" ./square-duct/"$dir"/0/"${f}"Dns

  done

  rm -rv ./square-duct/"$dir"/0/U
  rm -rv ./square-duct/"$dir"/0/R

done

for dir in ${ph_dirs};
do

  cp -rv ./periodic-hills/"$dir"/0/Udns ./periodic-hills/"$dir"/0/U
  cp -rv ./periodic-hills/"$dir"/0/Rdns ./periodic-hills/"$dir"/0/R
  mv -v ./periodic-hills/"$dir"/0/nut ./periodic-hills/"$dir"/0/nutRans

  calculateNut -case ./periodic-hills/"$dir"

  calculateRFV -case ./periodic-hills/"$dir"

  calculateRFVperp -case ./periodic-hills/"$dir"

  calculateRperp -case ./periodic-hills/"$dir"

  for f in ${fields};
  do

    mv -v ./periodic-hills/"$dir"/0/"$f" ./periodic-hills/"$dir"/0/"${f}"Dns

  done

  rm -rv ./periodic-hills/"$dir"/0/U
  rm -rv ./periodic-hills/"$dir"/0/R

done

