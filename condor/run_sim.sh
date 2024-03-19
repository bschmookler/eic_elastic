#!/bin/bash

DIR=$(cd -- $(dirname -- "${BASH_SOURCE[0]}") &> /dev/null  && pwd)
export LD_LIBRARY_PATH="${DIR}/install/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

echo "In directory $PWD!"

#Set Number of events
NEVENTS=100

#Input HepMC file
START=$((${NEVENTS}*$1))
FINISH=$((${NEVENTS}*($1+1)-1))
ln -s /gpfs02/eic/baraks/eic_elastic/hepmc/output/elas_gen_5_41_beameffects.hepmc input.hepmc

echo "Listing files before running:"
ls
echo ""

echo "Skipping first ${START} events!"
echo "Running ${NEVENTS} events!"
echo ""

source /opt/detector/setup.sh

#Run DIS events through npsim
npsim --compactFile $DETECTOR_PATH/epic_craterlake.xml --numberOfEvents ${NEVENTS} --skipNEvents ${START} --inputFiles input.hepmc  --outputFile output.edm4hep.root

#Run reconstruction
eicrecon -Ppodio:output_file=eicrecon_out.root -Pjana:nevents=${NEVENTS} -Pdd4hep:xml_files=epic_craterlake.xml output.edm4hep.root

#Clean up
mv output.edm4hep.root output_${START}_to_${FINISH}.edm4hep.root
mv eicrecon_out.root eicrecon_${START}_to_${FINISH}.root
unlink input.hepmc
rm tracking_geometry.obj

echo "Listing files after running:"
ls
echo ""

echo "Done!"

