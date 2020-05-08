echo "timepoint is "$1

matrix="methyl_"$1
x=$(seq -f "%g" -2 0.25 0)
echo $x

for i in $x 
do
	./ldak5.linux --sp $matrix \
	 --calc-kins-direct sens/${matrix}_${i}_alpha \
	 --SNP-data NO \
	 --ignore-weights YES \
	 --power ${i} \
	 --hwe-stand NO
done
