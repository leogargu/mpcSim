./test 

result=$(diff testSAM1_SAM_averaged.dat testSAM1_solution.dat -b -q | wc -l)

if [ "$result" == "0" ]; then
	echo "Test SAM 1 PASS"
else
	echo "Test SAM 1 FAIL"
fi

result=$(diff testCAM1_CAM_averaged.dat testCAM1_solution.dat -b -q | wc -l)
if [ "$result" == "0" ]; then
	echo "Test CAM 1 PASS"
else
	echo "Test CAM 1 FAIL"
fi