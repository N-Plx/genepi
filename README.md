# genepi

GENEPI Monte Carlo event generator. OSG version.

Describing electroproduction of photons and mesons off free and bound nucleons. Processes described : Bethe-Heitler, deeply virtual Compton scattering, pi0/eta production (and their decay to photons).

Documentation : `docs/` (physics model, sampling, options, output format, validation) and the original CLAS note https://misportal.jlab.org/ul/Physics/Hall-B/clas/viewFile.cfm/2009-024.pdf?documentId=554

Notable features :

- Events are selected with a keep/reject method. This means the output distributions already take into account the cross-sections.
- NH3/ND3 targets : events are randomly chosen to be on N or H/D. The output for N has a particle with ID 12. This is a trick, a fake ID to make it go through the simulation without interacting but still be able to recognize N events in the simulation output (gemc cannot track a generic C13/N13 ion from a LUND file: the real PDG code makes it exit).
- The vertex position and raster radius are chosen with --vx/--vy/--vz and --raster_x/--raster_y. The events are distributed uniformly inside an ellipse of centre (vx, vy) and half-axes (raster_x, raster_y), at z = vz.
- Events outside the physical limits (y, W2, nu, t, ...) are rejected, so the effective ranges can be tighter than the ones given on the command line, especially with Fermi motion in nuclei.

How to compile : needs ROOT (root-config, taken from $ROOTSYS if set, otherwise from the PATH), gfortran and make. On ifarm, for example:
> module use /scigroup/cvmfs/hallb/clas12/sw/modulefiles
>
> module load clas12
>
> make

The executable `genepi` is created in the top directory, next to gpd_table.dat. Check a build with:
> make -C tests && tests/run_smoke.sh

Clean :
> make clean

Run genepi :

> genepi [options]

    option   default  comment
    --trig     10000   number of output events
    --seed     auto    seed for random numbers, integer in [1, 2147483647]. If unspecified, a seed is generated
                       (random device, clock, process id, host name) and printed in the log. Use a unique
                       seed per job (e.g. from the grid job id) for reproducible productions.
    --docker   0       is ignored
    --Ebeam    10.6    incident e- momentum [GeV]
    --process  0       0 for DVCS, 1 for meson production.
    --meson    0       0 for pi0, 1 for eta
    --targ_A   1       Number of nucleons in the target
    --targ_Z   1       Number of protons in the target
    --x_min    0.001   Minimum x bjorken
    --x_max    1       Maximal x bjorken
    --y_min    0.001   Minimum beam energy fraction y=P.q/P.k  
    --y_max    1       Maximum beam energy fraction y=P.q/P.k    
    --Q2_min   1       Minimum photon virtuality, squared [GeV^2] Q2=−(k−k')^2
    --Q2_max   20      Maximum photon virtuality, squared [GeV^2] Q2=−(k−k')^2
    --W2_min   4       Minimum photon-nucleon invariant mass [GeV^2] W2= (P+q)^2
    --W2_max   50      Maximum photon-nucleon invariant mass [GeV^2] W2= (P+q)^2
    --nu_min   0.3     Minimum virtual photon energy [GeV] nu=Eb−E'
    --nu_max   11      Maximum virtual photon energy [GeV] nu=Eb−E'
    --t_min    -1.2    Minimum squared momentum transfer [GeV^2] t=(P−P')^2
    --t_max    0       Maximum squared momentum transfer [GeV^2] t=(P−P')^2
    --ycol_min -999    Minimum collider y ycol = (Q2+t)/(Q2+x*t)
    --ycol_max 0.025   Maximum collider y ycol = (Q2+t)/(Q2+x*t) 
    --vx       0.      Vertex position x
    --vy       0.      Vertex position y
    --vz       -3.     Vertex position z
    --raster_x 0.      Raster radius in the x direction [cm]
    --raster_y 0.      Raster radius in the y direction [cm]
    --output   genepi.dat  Output file (LUND format)
    --gamma_beam_min 1 DVCS only: minimum angle [deg] between the real photon and the beam
    --gamma_e_min    1 DVCS only: minimum angle [deg] between the real photon and the scattered electron
    --max_trials 0     Stop after this many trials (0 = no limit); the events kept so far are written.
                       Also replaces the default limit of 10,000,000 trials without any accepted event.
    
   where 
   - P,P' are the in/out nucleon 4-momenta
   - q the virtual photon 4-momentum
   - k/k' the in/out electron 4-momenta and Ebeam/E' their energies.

Invalid or inconsistent options (unknown option, non-numeric value, min >= max, ...) stop the program
with an error message. If no event is accepted after 10,000,000 trials (or --max_trials, if given), the
requested phase space is considered empty and the program stops (exit status 2).

The GPD table gpd_table.dat is read from $GENEPI if this variable is set, otherwise from the directory
of the executable, then from the current directory.

Exit status:

    0  success
    1  invalid input, missing GPD table, cannot open the output, or fatal JETSET error
    2  no event accepted (empty phase space)
    3  --max_trials reached after some, but fewer than --trig, events were kept
    4  error while writing the output file (e.g. disk full)
   
