#Clean previous output
rm -f ./SAM1/testSAM1_SAM_averaged.dat
rm -f ./SAM2/testSAM2_SAM_averaged.dat
rm -f ./SAM3/testSAM3_SAM_averaged.dat
rm -f ./SAM4/testSAM4_SAM_averaged.dat
rm -f ./CAM1/testCAM1_CAM_averaged.dat
rm -f ./CAM2/testCAM2_CAM_averaged.dat
rm -f ./CAMtoSAM/SAMconverted_testCAM1.dat
rm -f ./CAMtoSAM/SAMconverted_testCAM2.dat

#Run test
./test 


#Analyze results from test
result=$(diff ./SAM1/testSAM1_SAM_averaged.dat ./SAM1/testSAM1_solution.dat -b -q | wc -l)

if [ "$result" == "0" ]; then
	echo "Test SAM 1 ... PASS"
else
	echo "Test SAM 1 ... ***FAIL***"
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
	#printf "%s\n" "PASS" 
fi


result=$(diff ./SAM4/testSAM4_SAM_averaged.dat ./SAM4/testSAM4_solution.dat -b -q | wc -l)

if [ "$result" == "0" ]; then
	echo "Test SAM 4 ... PASS"
else
	echo "Test SAM 4 ... *FAIL*"
fi


result=$(diff ./CAM1/testCAM1_CAM_averaged.dat ./CAM1/testCAM1_solution.dat -b -q | wc -l)
if [ "$result" == "0" ]; then
	echo "Test CAM 1 ... PASS"
else
	echo "Test CAM 1 ... *FAIL*"
fi

result=$(diff ./CAM2/testCAM2_CAM_averaged.dat ./CAM2/testCAM2_solution.dat -b -q | wc -l)
if [ "$result" == "0" ]; then
	echo "Test CAM 2 ... PASS"
else
	echo "Test CAM 2 ... *FAIL*"
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