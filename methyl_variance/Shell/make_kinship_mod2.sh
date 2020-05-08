echo "timepoint is "$1

matrix="methyl_"$1

# Join weights
./ldak5.linux --sp $matrix \
 --join-weights mod2/$matrix \
 --SNP-data NO

./ldak5.linux --sp $matrix \
 --calc-kins-direct mod2/$matrix \
 --SNP-data NO \
 --weights mod2/$matrix/weights.short \
 --power -0.25 \
 --hwe-stand NO
