reset

set term postscript enhanced color
set output "~/Desktop/fE.ps"

path = "/Users/ashkandehghan/Desktop/SCFT_CODES/Graphene_Optimized/Graphene_Grand_Canonical_Optimized/RESULTS/fB2_0p8/FreeEnergyCalculations/"

set autoscale
set pointsize 2.0

set key top center font "Times,20" spacing 1.7

set xlabel '{/Symbol m}_H'  font "Times,32"
set ylabel 'L_x , L_y , L_z'  font "Times,32"

set xr [-14.0 : -12.7]

plot path."CAC.dat" using 2:($9/2.64) w lp lw 3 pt 4 title "L_x",\
path."CAC.dat" using 2:($10/4.56) w lp lw 3 pt 6 title "L_y",\
path."CAC.dat" using 2:($11/2.2) w lp lw 3 pt 8 title "L_z"

pause(-1)
