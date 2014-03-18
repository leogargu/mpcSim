README FILE FOR the current version(commit) of mpcSim

This is a file directory with the files for mpc Sim along with some examples (DATA). Previous tests and data are in mpcSim-1.0 in godzilla.

Further production DATA will be added in the future. Create a different repo with only production data?


Current priorities:
1.- Calculate fluctuations in momentum and kinetic energy, and export the data for plotting. fluctuations should reduce as approaching thermodynamics limit (n_part->inf)
3.- Calculate pressure
4.- Clean the code
5.- Automate equilibration?
6.- Implement radial averaging in the velocity profile
7.- Compare profile with NS results 
9.- For the velocity profile, calculate the slip predicted with teh SAM and the CAM method. Is the slip reduced with CAM? Compare to Bedkiha paper and MSc thesis.

Complete To-do list
2.- Find out where is the parameter lambda in Whitmers paper
3.- Modify (?) the streaming to check for particles scaping before the first out-of-boundary check (this is unnecessary if particles do not experience transversal forces)
	5.- Calculate total energy and momentum (and temperature?), and plot/export the data
6.- Rescale velocities to keep energy constant?
9.- Implement separate file with initial data to be read at the start of the simulation	
10.- Decouple the shift and encage procedure from collide. Call encage with shift first, then pass vector c_p (updated) to collide instead of shift	
	
Recently done:
1.- Use Git and stablish a workflow
2.- Use gprof to profile the code. Note: to show performance of all functions, prevent indentation with -fno-indent in the Makefile, Also, use -pg both compiling and linking
3.- stream-collide debugged using a previous version stored in Git.
7.- Calculate velocity profile. Export data. Visualise

Rejected:
8.- Rewrite the code so that the shift of the particles is actually done? (profiling shows this is not significant for performance)
1.- Switch to linear streaming instead of quartic solver (profiling 21/02/14 showed this is not a priority)
8.- Fix calculate_CAMdata_mv2 in mpc.c : SAM for temperature cannot be calculated here as the mean velocity over particles and across samples is needed:average.c can do it if mpc.c is modified to export the 3D velocity of all the particles


SUMMARY OF RECEND FINDINGS AND NEW TODOS

The standard deviations for equilibration (kinetic energy and momentum) cannot (probabably) be calculated on the fly. better run once, keep track of the seed, then redo
after having seen the plot. For this pupose, we need to export the K and J dat at each timestep, plot it, then find at what timestep it happens. Then rerun with same seed, and start exporting production data after that point.
During the production data, monitoring the K and J is not needed.

Averages must be done according to CAM to avoid bias. SAM calculation is still useful in two cases:
1.- for some quantities like rho and J, SAM=CAM  (prove this), all the time. This is useful when rho and/or J are used in hydrodynamic quantities. (hydrodynamic quantities are a function of the mechanical properties [Garcia2006])
2.- both averages are the same in equilibrium. This can be used to double-check that equilibrium has been attained.

Priority now is to export this equilibration data and automate its plotting. Stablish a switch in the call to the executable possibly to give it the seed, and to stablish production or test run.
Then, implement routines to export data for rho, J and K to feed it into the CAM routine (average).
Then, perform the appropriate routines to calculate derived hydrodynamic quantities such as temperature, from the above averages.
In this, possibly include the calculation of teh pressure (ideal gas).

Besides this, monitor the momentum absorbed by the wall at each timestep. Include somehow a way to see this in a thick slice, not only over the whole cylinder?
This will give an idea of the pressure on the walls.

Pending points are checking that teh temperature dip found by Garcia in the NATO report was obtained under the use of a thermostat. Compare with results from mpcSim.
Anotehr pending point is teh average of the velocity over the points of the same radius for a slice, and the comparison with the NS parabolic profile (or paraboloid).

Lastly but very important, estimate the slip when the fluid velociy is obtained with teh CAM average.

