#!/bin/bash
#SBATCH --job-name=trimGalore_JH		                                    #Job name
#SBATCH --partition=batch		                                        #Partition (queue) name
#SBATCH --ntasks=1			                                            #Single task job
#SBATCH --cpus-per-task=8	                                          #Number of cores per task
#SBATCH --mem=24gb			                                            #Total memory for job
#SBATCH --time=48:00:00  		                                        #Time limit hrs:min:sec
#SBATCH --output=/scratch/ahw22099/JH_GyneWorker/JH_trimmed_fq/log.%j			    #Standard output
#SBATCH --error=/scratch/ahw22099/JH_GyneWorker/JH_trimmed_fq/err.%j			    #Standard error log
#SBATCH --mail-user=ahw22099@uga.edu                                #Where to send mail -
#SBATCH --mail-type=END,FAIL                                        #Mail events (BEGIN, END, FAIL, ALL)

module load Trim_Galore/0.6.5-GCCcore-8.3.0-Java-11-Python-3.7.4

JH_trimmed_fq="/scratch/ahw22099/JH_GyneWorker/JH_trimmed_fq"
if [ ! -d $JH_trimmed_fq ]
then
mkdir -p $JH_trimmed_fq
fi

for file in /scratch/ahw22099/JH_GyneWorker/JH_raw_fq/*.gz
do
trim_galore --fastqc "$file" -o "$JH_trimmed_fq"
done
