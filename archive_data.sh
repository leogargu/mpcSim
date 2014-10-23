#echo $1
if [ "$1" == "" ]; then 
	echo "Error. Call as: bash archive_data.sh \"data title\" "
	echo "Aborting... "
	exit 0
fi

#Get latest git commit id
GIT_ID=$(git log --format="%H" -n 1)
#echo ${GIT_ID}

#Make data subfolder with current date and time
directory_name=$(date +"%Y.%m.%d_%H.%M")
mkdir data/$directory_name

#Move the data, pbs files and images from /experiments/ to the new folder
#mv ./experiments/* ./data/$directory_name

#write entry in /data/README.txt file with the name of the subfolder and a short title 
#(comments can be added later by hand) 
file="./data/README.txt"

#write git commit id into the file. [Note that line 5 must exist in the initial file for this to work!]
sed -i.bak "5,1i$directory_name - $1 \ngit commit id: $GIT_ID \ncomments: \n" $file

#print confirmation
echo "Completed. Add comments on this archived data by editing $file"
