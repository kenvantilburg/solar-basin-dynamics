driver4.jl is the main script used to integrate the solar system forwards 
and backwards with test particle injections.  driver4_XXX.jl are files 
used to increase the resolution of certain parts of the simulation, examine
how quickly particle orbits are lifted out of the sun, investigate GR effects,
etc.  Files like earth_hit_parr_XXX are used to determine how many particles 
are on sun/earth crossing orbits, and depend on external files like Output2.jl
or SolarBasinBasic.jl for formatting or other functions.  run_earth.sh is an
example SLURM script that calls driver4 to run a backwards simulation.