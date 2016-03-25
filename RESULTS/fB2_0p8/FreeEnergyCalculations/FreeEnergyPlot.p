reset

set term postscript enhanced color
set output "~/Desktop/fE.ps"



path = "/Users/ashkandehghan/Desktop/SCFT_CODES/Graphene_Optimized/Graphene_Grand_Canonical_Optimized/RESULTS/fB2_0p8/FreeEnergyCalculations/"

set autoscale
set pointsize 1.5

set key top right font " ,20" spacing 1.7

set xlabel '{/Symbol m}_H'  font " ,32"
set ylabel 'fE'  font " ,32"

set xr [-14.0 : -12.5]

plot path."Hom_Phase_2.dat" using 1:3 w lp lw 4 title "Homogenous",\
path."alphaBN.dat" using 2:7 w lp lw 4 title "{/Symbol a}-BN",\
path."Bilayer.dat" using 2:7 w lp lw 4 title "Single Layer {/Symbol a}-BN"

pause(-1)
