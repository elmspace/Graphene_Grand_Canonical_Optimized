#Script to cleanup data
clear

echo “Cleanup Script Has Started”

cd ./OMEGA/RUN_TIME_DATA
echo “Cleaning runtime omega files in ./OMEGA/RUN_TIME_DATA”
rm -r *
cd ..
cd ..

cd ./PHI
echo “Cleaning runtime phi files in ./PHI”
rm *
cd ..

cd ./MATLAB
echo “Cleaning density files in ./MATLAB”
rm -r *.dat
cd ..

cd ./RESULTS
echo “Cleaning fields in ./RESULTS”
rm -r *.dat
cd ..

echo “Cleanup Script Has Ended”
