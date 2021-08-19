#!/bin/bash
#SBATCH --job-name=JH_vcf		                                    #Job name
#SBATCH --partition=batch		                                        #Partition (queue) name
#SBATCH --ntasks=1			                                            #Single task job
#SBATCH --cpus-per-task=6                                          #Number of cores per task
#SBATCH --mem=24gb			                                            #Total memory for job
#SBATCH --time=48:00:00  		                                        #Time limit hrs:min:sec
#SBATCH --output=/scratch/ahw22099/JH_GyneWorker/STAR/Aligned_bam/log.%j			    #Standard output
#SBATCH --error=/scratch/ahw22099/JH_GyneWorker/STAR/Aligned_bam/err.%j			    #Standard error log
#SBATCH --mail-user=ahw22099@uga.edu                                #Where to send mail -
#SBATCH --mail-type=END,FAIL                                        #Mail events (BEGIN, END, FAIL, ALL)


SINV_GENOME="/scratch/ahw22099/JH_GyneWorker/UNIL_Sinv_3.0"
OUTDIR="/scratch/ahw22099/JH_GyneWorker/STAR/Aligned_bam"

samtools sort --threads 6  -o $OUTDIR/SRR7209532_trimmed.2pass.sorted.bam $OUTDIR/Aligned_bam/SRR7209532_trimmed.fq.gz.2pass.Aligned.sortedByCoord.out.bam
samtools index $OUTDIR/SRR7209532_trimmed.2pass.sorted.bam

SINV_GENOME="/scratch/ahw22099/JH_GyneWorker/UNIL_Sinv_3.0"
OUTDIR="/scratch/ahw22099/JH_GyneWorker/STAR/Aligned_bam"
module load BWA/0.7.17-GCC-8.3.0
bwa index $SINV_GENOME/GCF_016802725.1_UNIL_Sinv_3.0_genomic.fna

module load BCFtools/1.10.2-GCC-8.3.0

bcftools mpileup -O b -o $OUTDIR/SRR7209532_raw.bcf -f $SINV_GENOME/GCF_016802725.1_UNIL_Sinv_3.0_genomic.fna \
$OUTDIR/SRR7209532_trimmed.fq.gz.2pass.sorted.bam

bcftools call --threads 6 --ploidy 2 -m -v -o $OUTDIR/SRR7209532_trimmed.2pass.sorted.mpileup.vcf \
$OUTDIR/SRR7209532_raw.bcf

vcfutils.pl varFilter $OUTDIR/SRR7209532_trimmed.2pass.sorted.mpileup.vcf  > $OUTDIR/SRR7209532_trimmed.2pass.sorted.mpileup.final.vcf

bgzip $OUTDIR/SRR7209532_trimmed.2pass.sorted.mpileup.final.vcf

bcftools index $OUTDIR/SRR7209532_trimmed.2pass.sorted.mpileup.final.vcf





> $OUTDIR/SRR7209532_trimmed.2pass.sorted.mpileup.call.vcf
bcftools filter -Oz -e '%QUAL<40 || DP<10' $OUTDIR/SRR7209532_trimmed.2pass.sorted.mpileup.call.vcf > $OUTDIR/SRR7209532_trimmed.2pass.sorted.mpileup.call.filter.vcf
# Generate an IGV readable index file for your *vcf.gz file using `bcftools index` (BCFtools/1.10.2-GCC-8.3.0)
bcftools index $OUTDIR/SRR7209532_trimmed.2pass.sorted.mpileup.call.filter.vcf


cd $OUTDIR
for file in ./*.bam
  do
    samtools sort --threads 6 /scratch/ahw22099/JH_GyneWorker/STAR/Aligned_bam/SRR7209532_trimmed.fq.gz.2pass.Aligned.sortedByCoord.out.bam -o ./SRR7209532_trimmed.fq.gz.2pass.sorted.bam
    samtools index ./SRR7209532_trimmed.fq.gz.2pass.sorted.bam
  done

# Call variants with a (i) quality score of greater than 40, (ii) supported by more than 10 reads, (iii) with mapping quality greater than 60 for the the E. coli C600 genome using `bcftools mpileup`, `bcftools call`, and `bcftools filter` (BCFtools/1.10.2-GCC-8.3.0):
  # NOTE: Be sure to use multi-threading and match the `--threads` option in `bcftools mpileup` and and `bcftools call` to the `#SBATCH --cpus-per-task` SLURM value.
  # NOTE: Be sure to store your VCF file in *vcf.gz format and use the proper `--ploidy` value for a haploid bacteria when calling variants.
SINV_GENOME="/scratch/ahw22099/JH_GyneWorker/UNIL_Sinv_3.0"

# Construct a BWA index for the E. coli MG1655 refseq reference genome
module load BWA/0.7.17-GCC-8.3.0
bwa index -p Sinv_016802725.1_UNIL_Sinv_3.0 $SINV_GENOME/GCF_016802725.1_UNIL_Sinv_3.0
#
module load BCFtools/1.10.2-GCC-8.3.0
cd $OUTDIR
# for file in $OUTDIR/*.sorted.bam
#   do
bcftools mpileup -Oz --threads 6 --min-MQ 60 -f  ./SRR7209532_trimmed.2pass.sorted.bam > ./SRR7209532_trimmed.2pass.sorted.mpileup.vcf.gz
bcftools call -Oz -mv --threads 6 --ploidy 2 ./SRR7209532_trimmed.2pass.sorted.mpileup.vcf.gz > ./SRR7209532_trimmed.2pass.sorted.mpileup.call.vcf.gz
bcftools filter -Oz -e '%QUAL<40 || DP<10' ./SRR7209532_trimmed.2pass.sorted.mpileup.call.vcf.gz > ./SRR7209532_trimmed.2pass.sorted.mpileup.call.filter.vcf.gz
# Generate an IGV readable index file for your *vcf.gz file using `bcftools index` (BCFtools/1.10.2-GCC-8.3.0)
bcftools index ./SRR7209532_trimmed.2pass.sorted.mpileup.call.filter.vcf.gz
  done
