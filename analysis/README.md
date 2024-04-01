Analysis results:

1. [https://arxiv.org/abs/2207.04378](https://arxiv.org/abs/2207.04378): Results using the ROOT FOAM generator and fast detector simulations.
2. [Q2_range.py](analysis/Q2_range.py): For elastic scattering with a radiated photon, this code calculates and plots the true (i.e. vertex level) $Q^{2}$ as a function of the $x$ and $Q^{2}$ values based on measuring the scattered electron. The calculation used in the code come from equations 4-7 in [this paper](https://arxiv.org/abs/hep-ph/9906408v1). Note that for elastic scattering with a radiated photon, the true (i.e. vertex level) $x$ will be equal to 1.
3. In progress: Generator-level studies of DJANGOH elastic events including QED radiative effects.
4. [Elastic_reco.C](analysis/Elastic_reco.C): Analysis of ROOT FOAM events including beam-effects afterburner passed through the <i>ePIC</i> detector simulation and reconstruction. The resultant plots compare the generated and reconstructed kinematatics for the scattered electron and outgoing proton in both the lab frame and the (minimally-boosted) colinear beam reference frame. The reconstruction of $Q^{2}$ and $x$ is performed using several different methods.
