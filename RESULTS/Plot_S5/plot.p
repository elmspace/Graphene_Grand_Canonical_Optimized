

reset

#set term postscript enhanced color
#set output "~/Desktop/fE.ps"

set yr [-0.64 : -0.14]
set xr [0.05 : 0.21]

file = "/Users/ashkandehghan/Desktop/SCFT_CODES/Graphene_Optimized/Graphene_Grand_Canonical_Optimized/RESULTS/Plot_S5/fE.dat"
fileRef = "/Users/ashkandehghan/Desktop/SCFT_CODES/Graphene_Optimized/Graphene_Grand_Canonical_Optimized/RESULTS/Plot_S5/fE_Ref.dat"


set pointsize 2


plot file using 1:2 w lp pt 4 lw 3 title "{/Symbol a}-BN",\
fileRef using 1:3 w l lw 3 lc -1 title "{/Symbol a}-BN"

#file using 1:3 w lp pt 6 lw 3 title "C_{A/C}",\
#fileRef using 1:2 w l lw 3 lc -1 title "C_{A/C}",\

pause(-1)
