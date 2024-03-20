# eic_elastic
Generator and analysis for elastic electron-proton scattering at the EIC

Running ROOT-based generator
----------------------------
The ROOT-based generator in this repository -- [elas_gen_hepmc.C](elas_gen_hepmc.C) -- will generate elastic e-p events at EIC energies and write those events to a HepMC3 file. The events are generated using the Kelly Form-Factor parameterization, and the events are generated at the BORN level (no QED) radiation.

The generator can be run by doing:
```
root -l -b -q elas_gen_hepmc.C
```
The code will prompt the user to select the EIC beam energy combination, the upper and lower Q<sup>2</sup> limits, and the number of events to generate. The user will also be asked to select a generation option. There are two options:

1. <i>Uniform</i> generation. This will generate the scattered electron uniformly in solid angle in the proton beam rest frame. For each event, a cross-section weight is calculated and save to the 'weight' attribute in the HepMC file.
2. <i>FOAM</i> generation. This will generate the e-p events according the generator's cross-section model. In this case, no event-by-event weight needs to be applied. (In the output HepMC file, the 'weight' attribute is set to 1 for all events.)

The <i>FOAM</i> generation option is preferred for most studies.

Running DJANGOH generator
-------------------------
DJANGOH version 4.6.21 can generate elastic e-p events including QED radiation. (Note that generating these radiative elastic events will not work correctly in earlier versions of Djangoh.) 

To run DJANGOH v.4.6.21, clone and follow the instructions in [this external repository](https://github.com/tch285/DJANGOH). In this repository, the output format for the generated events has been modified so that the text output can be converted into the eic-smear(https://eic.github.io/software/eicsmear.html) and HepMC3 formats.

While you can compile DJANGOH locally, I compile and run using the [eic-shell environment](https://eic.github.io/tutorial-setting-up-environment/02-eic-shell/index.html), since this provides the pre-installed LHAPDF libraries and grids as well as the libraries needed for the eic-smear and HepMC3 file conversion.

Although the LHAPDF grids are not used for the elastic e-p running, you still need to set the location of the grids. When working in the <i>eic-shell</i>, you can link to the existing LHAPDF6 grids:
```
/cvmfs/sft.cern.ch/lcg/external/lhapdfsets/current
```

To run an elastic simulation with DJANGOH, use the [run_ep_rad_elas_5_41.sh](run_ep_rad_elas_5_41.sh) file in this local repository and do the following:
```
./run_ep_rad_elas_5_41.sh
```
This script runs DJANGOH, using the [ep_elas_rad_5x41.in](ep_elas_rad_5x41.in) input card to define the simulation parameters. The generated events are saved to a text file; a LOG file is also created, which contains the information about the simulation. The script will then run the [make_tree.C](make_tree.C) code to save the events to a ROOT file in the <i>eic-smear</i> format. Lastly, the script will run the [make_hepmc.C](make_hepmc.C) code to save the events in the HepMC3 format, which can then be used as input for the detector simulation.

Afterburner for beam-crossing and beam-smearing effects
--------------------------------------------------------
To include the EIC beam effects to the events created with the above generators, use this [afterburner](https://github.com/eic/afterburner) utility.

If running on the [eic-shell environment](https://eic.github.io/tutorial-setting-up-environment/02-eic-shell/index.html), this utility is included by default. To run the converter, do the following (for example):
```
abconv -o elas_gen_5_41_beameffects elas_gen_5_41.hepmc
```
This will create a file called ```elas_gen_5_41_beameffects.hepmc```, which can then be used as an input to the <i>ePIC</i> detector simulation.

Running generated events through the ePIC simulation
----------------------------------------------------
To run the 100 generated events through the detector simulation, access the [eic-shell environment](https://eic.github.io/tutorial-setting-up-environment/02-eic-shell/index.html).

Then simply do the following:
```
source /opt/detector/setup.sh

#Run DIS events through npsim
ln -s elas_gen_5_41_beameffects.hepmc input.hepmc #Just to make npsim command easier to read
npsim --compactFile $DETECTOR_PATH/epic_craterlake.xml --numberOfEvents 100 --inputFiles input.hepmc --outputFile output.edm4hep.root

#Run reconstruction
eicrecon -Ppodio:output_file=eicrecon_out.root -Pjana:nevents=100 -Pdd4hep:xml_files=epic_craterlake.xml output.edm4hep.root
```
The output ROOT file from this simulation can then be analyzed as discussed below.

To run events on a HTCondor-based batch farm, use the code in the [condor](condor) subdirectory.

Analysis
--------
Example analyses can be found in the [analysis](analysis) subdirectory.

