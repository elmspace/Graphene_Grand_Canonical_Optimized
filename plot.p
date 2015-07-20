

reset

#set term postscript enhanced color
#set output "~/Desktop/fE.ps"


file = "/Users/ashkandehghan/Desktop/SCFT_CODES/Graphene_Grand_Canonical/RESULTS/MOD1.dat"


set key at 0.1,2.2 font ",18" spacing 2.2

set xlabel '{/Symbol m}_H'  font ",22"
set ylabel '{/Symbol D}f'  font ",22"

#set xr [-20.0 : -13.4]


set pointsize 2

plot file using 2:($7-$8) w lp title "L_z" lw 4 pt 4 lc -1,\
"/Users/ashkandehghan/Desktop/SCFT_CODES/Graphene_Grand_Canonical/CALCULATIONS/kappa_0p5/MOD1.dat" using 2:($7-$8) w lp title "L_z" lw 4 pt 4 lc 1


pause(-1)

reset

#set term postscript enhanced color
#set output "~/Desktop/fE.ps"


file = "/Users/ashkandehghan/Desktop/SCFT_CODES/Graphene_Grand_Canonical/RESULTS/MOD1.dat"


set key at 0.12,2.25 font ",18" spacing 2.2

set xlabel '{/Symbol f}_H'  font ",22"
set ylabel 'L_x  L_y  L_z'  font ",22"

set xr [0.0 : 0.617565]


set pointsize 2

plot file using 6:($11/4.77) w lp title "L_z" lw 4 pt 4 lc -1,\
file using 6:($10/4.17) w lp title "L_y" lw 4 pt 6 lc -1,\
file using 6:($9/2.36) w lp title "L_x" lw 4 pt 8 lc -1

pause(-1)
