./test 

result=$(diff ./SAM1/testSAM1_SAM_averaged.dat ./SAM1/testSAM1_solution.dat -b -q | wc -l)

if [ "$result" == "0" ]; then
	echo "Test SAM 1 ... PASS"
else
	echo "Test SAM 1 ... *FAIL*"
fi

result=$(diff ./CAM1/testCAM1_CAM_averaged.dat ./CAM1/testCAM1_solution.dat -b -q | wc -l)
if [ "$result" == "0" ]; then
	echo "Test CAM 1 ... PASS"
else
	echo "Test CAM 1 ... *FAIL*"
fi

result=$(diff ./CAMtoSAM/testCAM1_SAMconverted.dat ./CAMtoSAM/testCAMtoSAM1_solution.dat -b -q | wc -l)
if [ "$result" == "0" ]; then
	echo "Test CAM to SAM 1 ... PASS"
else
	echo "Test CAM to SAM 1 ... *FAIL*"
fi

result=$(diff ./CAMtoSAM/testCAM2_SAMconverted.dat ./CAMtoSAM/testCAMtoSAM2_solution.dat -b -q | wc -l)
if [ "$result" == "0" ]; then
	echo "Test CAM to SAM 2 ... PASS"
else
	echo "Test CAM to SAM 2 ... *FAIL*"
fi