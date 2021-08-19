#!/bin/bash
#SBATCH --job-name=STAR_JH_1p2p	                                    #Job name
#SBATCH --partition=batch		                                        #Partition (queue) name
#SBATCH --ntasks=1			                                            #Single task job
#SBATCH --cpus-per-task=8	                                          #Number of cores per task
#SBATCH --mem=24gb			                                            #Total memory for job
#SBATCH --time=48:00:00  		                                        #Time limit hrs:min:sec
#SBATCH --output=/scratch/ahw22099/JH_GyneWorker/STAR/2p_out/log.%j			    #Standard output
#SBATCH --error=/scratch/ahw22099/JH_GyneWorker/STAR/2p_out/err.%j			    #Standard error log
#SBATCH --mail-user=ahw22099@uga.edu                                #Where to send mail -
#SBATCH --mail-type=END,FAIL                                        #Mail events (BEGIN, END, FAIL, ALL)

JH_trimmed_fq="/scratch/ahw22099/JH_GyneWorker/JH_trimmed_fq"
if [ ! -d $OWSF_trimmed_fq ]
then
mkdir -p $OWSF_trimmed_fq
fi

STAR_genome_idx="/scratch/ahw22099/JH_GyneWorker/alignment_genome_idx"
if [ ! -d $STAR_genome_idx ]
then
mkdir -p $STAR_genome_idx
fi
# create genome dir prior to running (add genome files from bghlab project dir)
SINV_GENOME="/scratch/ahw22099/JH_GyneWorker/UNIL_Sinv_3.0"
################## STAR ##################

module load STAR/2.7.2b-GCC-8.3.0
#GENOME GENERATE
time STAR \
--runMode genomeGenerate \
--runThreadN 8 \
--genomeSAindexNbases 13 \
--genomeDir $STAR_genome_idx \
--genomeFastaFiles $SINV_GENOME/GCF_016802725.1_UNIL_Sinv_3.0_genomic.fna \
--sjdbGTFfile $SINV_GENOME/GCF_016802725.1_UNIL_Sinv_3.0_genomic.gtf

##1st pass
STAR="/scratch/ahw22099/JH_GyneWorker/STAR"
if [ ! -d $MASTER ]
then
mkdir -p $MASTER
fi

FirstPass="/scratch/ahw22099/JH_GyneWorker/STAR/1p_out"
if [ ! -d $FirstPass ]
then
mkdir -p $FirstPass
fi

cd /scratch/ahw22099/JH_GyneWorker/JH_trimmed_fq/
for file in ./*fq.gz
do
run_STAR_1p () {
time STAR \
--readFilesCommand zcat \
--runThreadN 8 \
--genomeDir $STAR_genome_idx \
--readFilesIn $file \
--outFilterType BySJout \
--outFilterMultimapNmax 20 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNoverLmax 0.05 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--genomeLoad NoSharedMemory \
--outSAMtype BAM SortedByCoordinate \
--outSAMstrandField intronMotif \
--outSAMattrIHstart 0 \
--outFileNamePrefix /scratch/ahw22099/JH_GyneWorker/STAR/1p_out/"$file".1pass. \
--limitBAMsortRAM 15000000000
}
run_STAR_1p $file
done

##2nd pass
SecondPass="/scratch/ahw22099/JH_GyneWorker/STAR/2p_out"
if [ ! -d $SecondPass ]
then
mkdir -p $SecondPass
fi

cd /scratch/ahw22099/JH_GyneWorker/JH_trimmed_fq/
for file in ./*.fq.gz
do
run_STAR_2p () {
time STAR \
--readFilesCommand \
zcat --runThreadN 8 \
--quantMode GeneCounts \
--genomeDir $STAR_genome_idx \
--readFilesIn $file \
--outFilterType BySJout \
--outFilterMultimapNmax 20 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNoverLmax 0.05 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--genomeLoad NoSharedMemory \
--outSAMtype BAM SortedByCoordinate \
--outSAMstrandField intronMotif \
--outSAMattrIHstart 0 \
--outFileNamePrefix /scratch/ahw22099/JH_GyneWorker/STAR/2p_out/"$file".2pass. \
--limitBAMsortRAM 15000000000 \
--sjdbFileChrStartEnd \
$FirstPass/SRR7209532_trimmed.fq.gz.1pass.SJ.out.tab \
$FirstPass/SRR7209533_trimmed.fq.gz.1pass.SJ.out.tab \
$FirstPass/SRR7209534_trimmed.fq.gz.1pass.SJ.out.tab \
$FirstPass/SRR7209535_trimmed.fq.gz.1pass.SJ.out.tab \
$FirstPass/SRR7209536_trimmed.fq.gz.1pass.SJ.out.tab \
$FirstPass/SRR7209537_trimmed.fq.gz.1pass.SJ.out.tab \
$FirstPass/SRR7209538_trimmed.fq.gz.1pass.SJ.out.tab \
$FirstPass/SRR7209539_trimmed.fq.gz.1pass.SJ.out.tab \
$FirstPass/SRR7209540_trimmed.fq.gz.1pass.SJ.out.tab \
$FirstPass/SRR7209541_trimmed.fq.gz.1pass.SJ.out.tab \
$FirstPass/SRR7209542_trimmed.fq.gz.1pass.SJ.out.tab \
$FirstPass/SRR7209543_trimmed.fq.gz.1pass.SJ.out.tab \
$FirstPass/SRR7209544_trimmed.fq.gz.1pass.SJ.out.tab \
$FirstPass/SRR7209545_trimmed.fq.gz.1pass.SJ.out.tab \
$FirstPass/SRR7209546_trimmed.fq.gz.1pass.SJ.out.tab \
$FirstPass/SRR7209547_trimmed.fq.gz.1pass.SJ.out.tab \
$FirstPass/SRR7209548_trimmed.fq.gz.1pass.SJ.out.tab \
$FirstPass/SRR7209549_trimmed.fq.gz.1pass.SJ.out.tab \
$FirstPass/SRR7209550_trimmed.fq.gz.1pass.SJ.out.tab \
$FirstPass/SRR7209551_trimmed.fq.gz.1pass.SJ.out.tab
}
run_STAR_2p $file
done

### FINISHED 6/15/2021
