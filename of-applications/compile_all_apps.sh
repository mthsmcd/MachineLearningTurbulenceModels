#!/bin/bash

directories="calculateNut calculateRFV calculateRFVperp calculateRperp"

for dir in ${directories};
do

  echo -e "Compiling $dir"
  wclean ./$dir
  wmake ./$dir
  echo -e "Finished compiling $dir \n\n"

done