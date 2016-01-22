

reset

#set term postscript enhanced color
#set output "~/Desktop/fE.ps"

set autoscale

file = "/Users/ashkandehghan/Desktop/SCFT_CODES/Graphene_Optimized/Graphene_Grand_Canonical_Optimized/RESULTS/FreeEnergyCalculations/Hom_Phase_1.dat"
file2 = "/Users/ashkandehghan/Desktop/SCFT_CODES/Graphene_Optimized/Graphene_Grand_Canonical_Optimized/RESULTS/FreeEnergyCalculations/Hom_Phase_2.dat"
file3 = "/Users/ashkandehghan/Desktop/SCFT_CODES/Graphene_Optimized/Graphene_Grand_Canonical_Optimized/RESULTS/FreeEnergyCalculations/alphaBN.dat"
file4 = "/Users/ashkandehghan/Desktop/SCFT_CODES/Graphene_Optimized/Graphene_Grand_Canonical_Optimized/RESULTS/FreeEnergyCalculations/Bilayer.dat"

plot file using 1:3 w lp title "fE",\
file2 using 1:3 w lp title "fE",\
file3 using 2:7 w lp title "fE",\
file4 using 2:7 w lp title "fE"


pause(-1)
