#!/bin/tcsh

setenv EIC_SHELL_PREFIX '/gpfs02/eic/baraks/epic/local'
setenv SINGULARITY_BINDPATH '/gpfs02,/gpfs01,/gpfs'
/usr/bin/singularity exec $EIC_SHELL_PREFIX/lib/jug_xl-nightly /bin/bash -c "$argv[1-]"
