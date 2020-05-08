echo "timepoint is "$1

#matrix="methyl_FOM"
#cut="FALSE"
#start_section=54
#end_section=423
matrix="methyl_"$1
cut=$2
start_section=$3
end_section=$4

if [ "$cut" == "TRUE" ]; then
	echo "Cutting weights"

	./ldak5.linux --sp $matrix \
	 --cut-weights mod2/$matrix \
	 --SNP-data NO

	echo "Finished cutting weights"
fi

section=$3

# Calculate weights
for i in `seq $start_section $end_section`;
do
	./ldak5.linux --sp $matrix \
	 --calc-weights mod2/$matrix \
	 --section $i \
	 --SNP-data NO 
done

echo "Finished calculating weights"
