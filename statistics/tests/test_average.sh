./test 

result=$(diff ./SAM1/testSAM1_SAM_averaged.dat ./SAM1/testSAM1_solution.dat -b -q | wc -l)

if [ "$result" == "0" ]; then
	echo "Test SAM 1 ... PASS"
else
	echo "Test SAM 1 ... *FAIL*"
fi

result=$(diff ./SAM2/testSAM2_SAM_averaged.dat ./SAM2/testSAM2_solution.dat -b -q | wc -l)

if [ "$result" == "0" ]; then
	echo "Test SAM 2 ... PASS"
else
	echo "Test SAM 2 ... *FAIL*"
fi

result=$(diff ./SAM3/testSAM3_SAM_averaged.dat ./SAM3/testSAM3_solution.dat -b -q | wc -l)

if [ "$result" == "0" ]; then
	echo "Test SAM 3 ... PASS"
else
	echo "Test SAM 3 ... *FAIL*"
fi

result=$(diff ./CAM1/testCAM1_CAM_averaged.dat ./CAM1/testCAM1_solution.dat -b -q | wc -l)
if [ "$result" == "0" ]; then
	echo "Test CAM 1 ... PASS"
else
	echo "Test CAM 1 ... *FAIL*"
fi

result=$(diff ./CAMtoSAM/SAMconverted_testCAM1.dat ./CAMtoSAM/testCAMtoSAM1_solution.dat -b -q | wc -l)
if [ "$result" == "0" ]; then
	echo "Test CAM to SAM 1 ... PASS"
else
	echo "Test CAM to SAM 1 ... *FAIL*"
fi

result=$(diff ./CAMtoSAM/SAMconverted_testCAM2.dat ./CAMtoSAM/testCAMtoSAM2_solution.dat -b -q | wc -l)
if [ "$result" == "0" ]; then
	echo "Test CAM to SAM 2 ... PASS"
else
	echo "Test CAM to SAM 2 ... *FAIL*"
fi