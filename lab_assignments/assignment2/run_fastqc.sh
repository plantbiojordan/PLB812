#!/bin/bash --login
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=128GB
#SBATCH --job-name fastqc_la2
#SBATCH --mail-user=manchego@msu.edu
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=ALL

# directory with fastq files
cd /mnt/home/manchego/PLB812/data/ncbi/public/sra

# bringing in dependencies and packages from conda env plb812
export PATH="${HOME}/miniconda3/envs/plb812/bin:${PATH}"
export LD_LIBRARY_PATH="${HOME}/miniconda3/envs/plb812/lib:${LD_LIBRARY_PATH}"

# list of fastq files in /sra
fqfiles = "SRR5448192.fastq SRR5448193.fastq SRR5448194.fastq SRR5448195.fastq \
SRR5448196.fastq SRR5448197.fastq SRR5448198.fastq SRR5448199.fastq \
SRR5448200.fastq SRR5448201.fastq SRR5448202.fastq SRR5448203.fastq"

# for loop to produce .html output files from fastqc
for thing in ${fqfiles}
do

fastqc ${thing}

done


# copying those files to PLB812/lab_assignments/assignment2
cd /mnt/home/manchego/PLB812

cp /mnt/home/manchego/PLB812/data/ncbi/public/sra/*html \
/mnt/home/manchego/PLB812/lab_assignments/assignment2/
