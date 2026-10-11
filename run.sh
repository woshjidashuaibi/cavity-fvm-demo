#!/bin/bash  
echo "Compiling cavity_fvm.c ..."  
gcc cavity_fvm.c -o cavity_fvm -lm  
echo "Run cavity FVM lid-driven cavity demo"  
./cavity_fvm 
