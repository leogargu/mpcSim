#Number of tests [update this when writing new tests]
#[remember also to update test_average.c and compile it]
numSAM=5
numCAM=4
numCAMtoSAM=2

#Auxiliary strings
filename_SAM_end="_SAM_averaged.dat"
filename_CAM_end="_CAM_averaged.dat"

#Clean previous output
#SAM tests
for(( i=1; i<=$numSAM; i++  )); do
	filename_start="./SAM$i/testSAM$i"
	filename_full=$filename_start$filename_SAM_end
	rm -f $filename_full
done
#CAM tests
for(( i=1; i<=$numCAM; i++  )); do
	filename_start="./CAM$i/testCAM$i"
	filename_full=$filename_start$filename_CAM_end
	rm -f $filename_full
done
#CAMtoSAM tests
for(( i=1; i<=$numCAMtoSAM; i++  )); do
	filename_full="./CAMtoSAM/SAMconverted_testCAM$i.dat"
	rm -f $filename_full
done


#Run test
./test 


#Analyze results from tests
#SAM
for(( i=1; i<=$numSAM; i++)); do
	directory="./SAM$i/testSAM$i"
	fileout_end="_SAM_averaged.dat"
	fileout=$directory$fileout_end
	filesol_end="_solution.dat"
	filesol=$directory$filesol_end
	
	result=$(diff "$fileout" "$filesol" -b -q | wc -l)

	if [ "$result" == "0" ]; then
		echo "Test SAM $i ... PASS"
	else
		echo "Test SAM $i ... *FAIL*"
	fi
done

#CAM
for(( i=1; i<=$numCAM; i++)); do
	directory="./CAM$i/testCAM$i"
	fileout_end="_CAM_averaged.dat"
	fileout=$directory$fileout_end
	filesol_end="_solution.dat"
	filesol=$directory$filesol_end
	
	result=$(diff "$fileout" "$filesol" -b -q | wc -l)

	if [ "$result" == "0" ]; then
		echo "Test CAM $i ... PASS"
	else
		echo "Test CAM $i ... *FAIL*"
	fi
done

#CAMtoSAM 
for(( i=1; i<=$numCAMtoSAM; i++)); do
	fileout="./CAMtoSAM/SAMconverted_testCAM$i.dat"
	ending="_solution.dat"
	sol="./CAMtoSAM/testCAMtoSAM$i"
	filesol=$sol$ending
	
	result=$(diff "$fileout" "$filesol" -b -q | wc -l)
	if [ "$result" == "0" ]; then
		echo "Test CAM to SAM $i ... PASS"
	else
		echo "Test CAM to SAM $i ... *FAIL*"
	fi
done



#NOTES
#result=$(diff <(tail -n +2 ./SAM1/testSAM1_SAM_averaged.dat) <(tail -n +2 ./SAM1/testSAM1_solution.dat) -b -q | wc -l)
