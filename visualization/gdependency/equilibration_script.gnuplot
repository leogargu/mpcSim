set terminal postscript
set output "equilibration_energy.ps"
set xlabel "timesteps"
set ylabel "E_k/particle"
set title "Kinetic Energy equilibration"
set key off
plot "equilibration.dat" using 1:5 with lines
set output "equilibration_temperature.ps"
set xlabel "timesteps"
set ylabel "T (whole system)"
set title "Temperature equilibration"
set key off
plot "equilibration.dat" using 1:6 with lines
set output "equilibration_momentum.ps"
set xlabel "timesteps"
set ylabel "J/particle"
set title "Momentum equilibration"
set key above box
show key
plot "equilibration.dat" using 1:2 with lines, "equilibration.dat" using 1:3 with lines lc rgb "red", "equilibration.dat" using 1:4 with lines lc rgb "blue"
