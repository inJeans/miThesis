#!/bin/bash
# Majorana problem flip
BX=1.0e-7  # transverse magnetic field (T)
C=2.5e-3  # rate of change of z magnetic field (T/s)
DT=1.0e-6  # simulation time step (s)
FILENAME="majorana_flip_data"  # output data filename
FIGURENAME="flip"  # output figure file prefix
FIGUREDIR="../../gfx/Chapter-4/"  # output figure directory

python majorana_problem.py -bx $BX -c $C -dt $DT -o $FILENAME
python majorana_data_crunch.py -i $FILENAME -o $FIGURENAME -d $FIGUREDIR