reset

file = "/Users/ashkandehghan/Desktop/SCFT_CODES/Graphene_Optimized/Graphene_Grand_Canonical_Optimized/PHI/phiyz.dat"

set pm3d
set iso 100
set samp 100
set palette model RGB
set dgrid3d 32,32,1
set size square
unset border
unset xtics
unset ytics
unset colorbox
set pm3d flush begin ftriangles scansforward interpolate 10,5
   
   unset key
   unset sur
   set hidden3d
   set view map 
   set autoscale
   set size square
   
   splot file using 1:2:4
