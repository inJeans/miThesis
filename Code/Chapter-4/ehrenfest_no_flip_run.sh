#!/bin/bash
# Ehrenfest method
BX=1.45e-6  # transverse magnetic field (T)
DBZ=2.5  # rate of change of z magnetic field (T/s)
DT=1.0e-7  # simulation time step (s)
FILENAME="ehrenfest_data"  # output data filename
FIGURENAME="ehrenfest"  # output figure file prefix
FIGUREDIR="../../gfx/Chapter-4/"  # output figure directory

python ehrenfest_method.py -bx $BX -dbz $DBZ -dt $DT -o $FILENAME
python ehrenfest_data_crunch.py -i $FILENAME -o $FIGURENAME -d $FIGUREDIR