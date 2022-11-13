# genepi

GENEPI Monte Carlo event generator.

Describing electroproduction of photons and mesons off free and bound nucleons. Processes described : Bethe-Heitler, deeply virtual Compton scattering, pi0/eta production (and their decay to photons).

Documentation : https://misportal.jlab.org/ul/Physics/Hall-B/clas/viewFile.cfm/2009-024.pdf?documentId=554

Additional features :

> Events are selected with a keep/reject method. This means the output distributions already take into account the cross-sections.
> NH3/ND3 targets : events are randomly chosen to be on N or H/D. The output for N has a particle with ID 12. This is a trick, a fake ID to make it go through the simulation without interacting but still be abl\
e to recognize N events in the simulation output.
> Phi meson use the pi0 cross-sections.


How to compile (on ifarm, on a clean environment) :
> mkdir bin //Create an empty bin directory where the executable will live
>
> source setenv.txt
> 
> make 

Clean executables :
> make clean

Run genepi :

> bin/genepi.exe input_file.input

An example of input_file can be found in 

> input_files/
