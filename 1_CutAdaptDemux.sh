#!/bin/bash -e
#SBATCH -A uoo02328
#SBATCH -J CutAdapt-emtDNA
#SBATCH --time 12:00:00
#SBATCH -N 1
#SBATCH -c 24
#SBATCH -n 1
#SBATCH --mem=48G
#SBATCH --partition=large
#SBATCH --mail-user=clare.adams@postgrad.otago.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --output CutAdapt-emtDNA.%j.out # CHANGE map1 part each run
#SBATCH --error CutAdapt-emtDNA.%j.err # CHANGE map1 part each run

module load cutadapt/2.10-gimkl-2020a-Python-3.8.2

cd /nesi/nobackup/uoo02328/clare/Paua-e-mtDNA/ || exit

RAW=/nesi/nobackup/uoo02328/clare/Paua-e-mtDNA/RAW/

echo "making CutAdapt folder and loading in CutAdapt"

mkdir CutAdapt

cd CutAdapt

module load FastQC/0.11.9
fastqc /nesi/nobackup/uoo02328/clare/Paua-e-mtDNA//RAW/*.gz

echo "Running cutAdapt demultiplexing script"

cutadapt \
    -e 0.15 --no-indels \
    -m 100 \
    -g file:${RAW}/CutAdapt-F-barcodes.fasta \
    -G file:${RAW}/CutAdapt-R-barcodes.fasta \
    -o trimmed-{name1}-{name2}.1.fastq.gz -p trimmed-{name1}-{name2}.2.fastq.gz \
    ${RAW}/6119-P1-0-1_S1_L001_R1_001.fastq.gz ${RAW}/6119-P1-0-1_S1_L001_R2_001.fastq.gz

echo "Done! Ready to download to computer for renaming"
