#!/usr/bin/bash

echo "-----------------------------------"
echo "Running DJANGOH Simulation for ep Collider!!!"
echo "..."
echo ""

OUTFILE1=outfiles/djangoh.rad_elas_5x41_evt.dat
if test -f "$OUTFILE1"; then
	rm -f "$OUTFILE1"
fi
OUTFILE2=outfiles/djangoh.rad_elas_5x41_out.dat
if test -f "$OUTFILE2"; then
	rm -f "$OUTFILE2"
fi
OUTFILE3=outfiles/djangoh.rad_elas_5x41_smp.dat
if test -f "$OUTFILE3"; then
        rm -f "$OUTFILE3"
fi

./djangoh < ep_elas_rad_5x41.in > logfiles/ep_elas_rad_5x41.log

echo "Completed Simulation!!!"
echo ""

echo "Making Output ROOT File..."
root -l -b -q 'make_tree.C("djangoh.rad_elas_5x41_evt.dat")'
echo "Done!!!"

echo "Making Output HEPMC File..."
root -l -b -q 'make_hepmc.C("djangoh.rad_elas_5x41_evt.root")'
echo "-----------------------------------"

echo "-----------------------------------"

