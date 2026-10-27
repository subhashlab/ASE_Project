#!/bin/bash

#### Generate job scripts
#./generate_batch_scripts.sh 1 2 5gb
#./generate_batch_scripts.sh 1 16 15gb
#### Submit multipe jobs
#find . -name '*.sbatch' | xargs -I% qsub %

### Sbatch parameters
numnode=$1
cpus=$2
#time=$3
cputype="cpu-test" #"cpu-parallel"
memreq=$3

SCRIPT_PATH="/iitjfs/home/2024rbb1003/projects/ASE_NB_project/RIPSeq/scripts";
SAMPLES=($(cat "$SCRIPT_PATH/samplenames_input.txt"))


for i in ${SAMPLES[@]}

do

	IFS=';' read -ra arr <<< "$i"; 

	echo "#!/bin/bash -l" >> $SCRIPT_PATH/bscripts/${arr[0]}.sbatch
	echo "#PBS -q $cputype" >> $SCRIPT_PATH/bscripts/${arr[0]}.sbatch
	echo "#PBS -l nodes=$numnode:ppn=$cpus" >> $SCRIPT_PATH/bscripts/${arr[0]}.sbatch
	echo "#PBS -N ${arr[0]}_${arr[1]}_${arr[2]}" >> $SCRIPT_PATH/bscripts/${arr[0]}.sbatch
	echo "#PBS -l mem=$memreq" >> bscripts/${arr[0]}.sbatch
	echo "#PBS -o ${arr[0]}_out.txt" >> bscripts/${arr[0]}.sbatch
	echo "#PBS -e ${arr[0]}_err.txt" >> bscripts/${arr[0]}.sbatch

	#echo "echo Removing Adapter Sequence ${arr[0]}:${arr[1]}:${arr[2]} \`date\`;" >> bscripts/${arr[0]}.sbatch;
	#		echo "$SCRIPT_PATH/run_cutadapt.sh ${arr[0]} ;" >> $SCRIPT_PATH/bscripts/${arr[0]}.sbatch;
	#		echo "" >> $SCRIPT_PATH/bscripts/${arr[0]}.sbatch;

	#echo "echo Adapter Sequence Successfully Removed ${arr[0]}:${arr[1]}:${arr[2]} \`date\`;" >> bscripts/${arr[0]}.sbatch;
	#echo "wait" >> bscripts/${arr[0]}.sbatch;
	
	echo "echo Aligning started ${arr[0]}:${arr[1]}:${arr[2]} \`date\`;" >> bscripts/${arr[0]}.sbatch;
	#### Aligning (If the fastq files are not cleaned, please include adapter cleaning step before aligning)

			echo "$SCRIPT_PATH/run_STAR_align_Input_TE.sh ${arr[0]} $cpus ;" >> $SCRIPT_PATH/bscripts/${arr[0]}.sbatch;
			echo "" >> $SCRIPT_PATH/bscripts/${arr[0]}.sbatch;

	echo "echo Alignment completed Successfully ${arr[0]}:${arr[1]}:${arr[2]} \`date\`;" >> bscripts/${arr[0]}.sbatch;
	echo "wait" >> bscripts/${arr[0]}.sbatch;
	

done
