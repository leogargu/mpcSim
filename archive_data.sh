#Checking calling parameter(s) are valid
if [ "$1" == "" ]; then 
	echo "Error. Call as: bash archive_data.sh \"data title\" "
	echo "Aborting... "
	exit 0
fi

#Get latest git commit id
GIT_ID=$(git log --format="%H" -n 1)
#echo ${GIT_ID}

#Make data subfolder with current date and time
#[Careful with spaces before/after the = sign!!!!]
target_folder=$(date +"%Y.%m.%d_%H.%M")
target_directory=./data/$target_folder
mkdir $target_directory

#Determine origin of files to be archived
origin_directory=./experiments/

#Move everything to the new folder
# First copy the files...
rsync -aP $origin_directory $target_directory
# (Run command above again if the transfer stalls r is interrupted: rsync will resume from teh point where it stoppped)

#sanity check
origin_num=$(find $origin_directory -name "*" | wc -l)
target_num=$(find $target_directory -name "*" | wc -l)

if [ $origin_num == $target_num ]; then
	echo "Successfully copied files from "$origin_directory" into "$target_directory
else
	echo "Error while copying files. Exiting..."
	exit 1
fi
#...Then delete files inorigin directory, if required. Note rm * can spit a "argument list too long" error, so...
read -p "DELETE files in "$origin_directory" (yes/no) ? " CONT
if [ "$CONT" == "yes" ]; then
	find $origin_directory -maxdepth 1 -name "*.*" -print0 | xargs -0 rm
else
	echo "Files in "$origin_directory" have not been removed"
fi

#http://www.pronego.com/helpdesk/knowledgebase.php?article=59

#write entry in /data/README.txt file with the name of the subfolder and a short title 
#(comments can be added later by hand) 
file="./data/README.txt"

#write git commit id into the file. [Note that line 5 must exist in the initial file for this to work!]
sed -i.bak "5,1i$target_folder - $1 \ngit commit id: $GIT_ID \ncomments: \n" $file

#print confirmation
echo "COMPLETED. Add comments on this archived data by editing $file"

###############################################################################
#NOTES
#if [ rm $origin_directory* ]; then
#		echo "Files removed from "$origin_directory
#	else
#		echo "Too many files in "$origin_directory", removing one by one:"
		
		#for i in $origin_directory*; do
	#		echo "Deleting : $i";
	#		rm $i;
	#	done
#	fi
