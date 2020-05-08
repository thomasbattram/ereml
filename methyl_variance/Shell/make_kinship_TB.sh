echo "timepoint is "$1

matrix="methyl_"$1

./ldak5.linux --sp $matrix \
 --calc-kins-direct $matrix \
 --SNP-data NO \
 --ignore-weights YES \
 --power -0.25 \
 --hwe-stand NO
