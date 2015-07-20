

reset

#set term postscript enhanced color
#set output "~/Desktop/fE.ps"


file_x = "/Users/ashkandehghan/Desktop/SCFT_CODES/Graphene_Grand_Canonical/PHI/phix.dat"


set pointsize 2

plot file_x using 1:2 w lp title "A" lw 4 pt 4 lc 1,\
file_x using 1:3 w lp title "C" lw 4 pt 4 lc 3,\
file_x using 1:($4+$5+$6) w lp title "B" lw 4 pt 4 lc 2,\
file_x using 1:7 w lp title "Hom" lw 4 pt 4 lc 4


