# genepi

GENEPI Monte Carlo event generator.

Describing electroproduction of photons and mesons off free and bound nucleons. Processes described : Bethe-Heitler, deeply virtual Compton scattering, pi0/eta production (and their decay to photons).

Documentation : https://misportal.jlab.org/ul/Physics/Hall-B/clas/viewFile.cfm/2009-024.pdf?documentId=554

Additional features :

- Events are selected with a keep/reject method. This means the output distributions already take into account the cross-sections.
- NH3/ND3 targets : events are randomly chosen to be on N or H/D. The output for N has a particle with ID 12. This is a trick, a fake ID to make it go through the simulation without interacting but still be able to recognize N events in the simulation output. 
- Phi meson use the pi0 cross-sections. 

How to compile (on ifarm, on a clean environment) :
> setenv ROOTSYS /apps/root/6.10.02/root
> 
> setenv ROOTLIB /apps/root/6.10.02/root/lib
> 
> setenv LD_LIBRARY_PATH <span>$</span>{PATH}:${ROOTSYS}/lib
>
> setenv PATH <span>$</span>{PATH}:${ROOTSYS}/bin
> 
> make 

Clean executables :
> make clean

Run genepi :

> genepi [options]

    option   default  comment
    --trig     10000   number of output events
    --seed     time    seed for random numbers. if unspecified, uses current time information
    --docker   0       is ignored
    --Ebeam    10.6    incident e- momentum [GeV]
    --process  0       0 for DVCS, 1 for meson production.
    --meson    0       0 for pi0, 1 for eta, 2 for phi
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
    
   where 
   - P,P' are the in/out nucleon 4-momenta
   - q the virtual photon 4-momentum
   - k/k' the in/out electron 4-momenta and Ebeam/E' their energies.
   
