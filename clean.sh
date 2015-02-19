#!/bin/sh

scons -c 

for file in $(find . -maxdepth 1 -type f -name "*.out" )
do
  rm -rf $file 
done

for file in $(find . -maxdepth 1 -type f -name "*.frag?" )
do
  rm -rf $file 
done

for file in $(find . -maxdepth 1 -type f -name "*.log" )
do
  rm -rf $file 
done
