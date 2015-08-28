

reset

#set term postscript enhanced color
#set output "~/Desktop/fE.ps"


file = "/Users/ashkandehghan/Desktop/SCFT_CODES/Graphene_Optimized/Graphene_Grand_Canonical_Optimized/PHI/phix.dat"



plot file using 1:2 w lp title "A",\
file using 1:3 w lp title "C",\
file using 1:4 w lp title "B1",\
file using 1:5 w lp title "B2",\
file using 1:6 w lp title "B3",\
file using 1:7 w lp title "B4"


