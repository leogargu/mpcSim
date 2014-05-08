ncells=15
maxocc=50
maxnum=2

echo "This test will erase .dat, .out files in the current directory"

read -p "Continue (y/n)?" CONT
if [ "$CONT" == "y" ]; then
	echo "Test starting...";
	for i in {0..3}
	do
		rm *.dat *.out
		./test $ncells $maxocc $maxnum $i > screen.out
		tail -n +2 *.dat > file.out
		echo "Checking test case" $i
		diff screen.out file.out -b
	done

	echo "Test complete."
else
  echo "Aborted.";
fi




