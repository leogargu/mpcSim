README FILE FOR the current version(commit) of mpcSim

This is a file directory with the files for mpc Sim along with some examples (DATA). Previous tests and data are in mpcSim-1.0 in godzilla.

Further production DATA will be added in the future. Create a different repo with only production data?


Current priorities:
1.- Calculate fluctuations in momentum and kinetic energy, and export the data for plotting
2.- Calculate temperature, check it is constant
3.- Calculate pressure
4.- Clean the code
5.- Automate equilibration?
6.- Implement radial averaging in the velocity profile
7.- Compare profile with NS results 
8.- Build a debug system for teh stream-collide laternating algorith, Initiated in Matehmatica

Complete To-do list

	1.- Switch to linear streaming instead of quartic solver (profiling 21/02/14 showed this is not a priority)
2.- Find out where is the parameter lambda in Whitmers paper
3.- Modify (?) the streaming to check for particles scaping before the first out-of-boundary check (this is unnecessary if particles do not experience transversal forces)
5.- Calculate total energy and momentum (and temperature?), and plot/export the data
6.- Rescale velocities to keep energy constant?
7.- Calculate velocity profile. Export data. Visualise
	8.- Rewrite the code so that the shift of the particles is actually done? (profiling shows this is not significant for performance)
9.- Implement separate file with initial data to be read at the start of the simulation	
10.- Decouple the shift and encage procedure from collide. Call encage with shift first, then pass vector c_p (updated) to collide instead of shift	
	
Recently done:
1.- Use Git and stablish a workflow
2.- Use gprof to profile the code. Note: to show performance of all functions, prevent indentation with -fno-indent in the Makefile, Also, use -pg both compiling and linking
