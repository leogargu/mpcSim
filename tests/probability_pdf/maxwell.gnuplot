set terminal postscript
set output "maxwell.ps"
set xlabel "velocity"
set ylabel "P"
set title "Maxwell-Boltzmann RNG"
set key off
f(x) = sqrt(1/(2*pi)) * exp(-x*x*0.5)
plot [x=-5:5] f(x)  with lines, "maxwell.dat" using 1:2 with boxes


