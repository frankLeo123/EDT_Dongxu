#!/bin/bash 
./main $1 $2 $3 $4
gprof -b -p main
