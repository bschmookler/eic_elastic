# eic_elastic
Generator and analysis for elastic electron-proton scattering at the EIC

Running ROOT-based generator
----------------------------
The ROOT-based generator in this repository will generate elastic e-p events at EIC energies and write those events to a HepMC3 file. The events are generated using the Kelly Form-Factor parameterization, and the events are generated at the BORN level (no QED) radiation.

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
In progress...

Afterburner for beam-crossing and beam-smearing effects
--------------------------------------------------------
To include the EIC beam effects to the events created with the above generators, use this [afterburner](https://github.com/eic/afterburner) utility.

Analysis
--------
In progress...

