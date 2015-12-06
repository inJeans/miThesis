#!/bin/bash
# Majorana problem no flip
BX=2.0e-7  # transverse magnetic field (T)
C=2.5e-3  # rate of change of z magnetic field (T/s)
DT=1.0e-6  # simulation time step (s)
FILENAME="majorana_no_flip_data"  # output data filename
FIGURENAME="no_flip"  # output figure file prefix
FIGUREDIR="../../gfx/Chapter-4/"  # output figure directory

python majorana_problem.py -bx $BX -c $C -dt $DT -o $FILENAME
python majorana_data_crunch.py -i $FILENAME -o $FIGURENAME -d $FIGUREDIR