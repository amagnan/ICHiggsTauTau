#!/bin/bash

for dir in cards_run2ana_170217/
  do cd $dir
  echo $dir:

  for subdir in enu munu taunu mumu ee qcd nunu
    do cd $subdir
    echo $subdir:

    mv vbfhinv_${subdir}_13TeV.txt tmp.txt
    sed "s/\-\/\-/\-/g" tmp.txt > tmp1.txt
    sed "s/\-\//1\//g"            tmp1.txt > tmp2.txt
    sed "s/\/\-/\/1/g"                       tmp2.txt > vbfhinv_${subdir}_13TeV.txt
    rm tmp*
    cd ../
    done

  done
