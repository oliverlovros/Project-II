#!/bin/bash
# Run the program.

read -p 'How many processors do you want? ' processor
time mpirun -n $processor program.exe