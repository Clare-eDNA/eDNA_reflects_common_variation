#!/bin/bash -e
#SBATCH -A uoo02328
#SBATCH -J maxEE
#SBATCH --time 1:00:00
#SBATCH -N 1
#SBATCH -c 8
#SBATCH -n 1
#SBATCH --mem=15G
#SBATCH --partition=large
#SBATCH --mail-user=clare.adams@postgrad.otago.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --output ee.%j.out # CHANGE map1 part each run
#SBATCH --error ee.%j.err # CHANGE map1 part each run


module load USEARCH/11.0.667-i86linux32


cd /nesi/nobackup/uoo02328/clare/Paua-e-mtDNA/CutAdapt4/Samples/PlaywithMaxEE/maxEE-05 || exit


#Filter trimmed sequences for a max expected error rate
for sample in *.309.fastq
do
	out=$(basename "$sample" .309.fastq)
	echo "$out"
	usearch -fastq_filter "$sample" -fastq_maxee 0.05 -relabel "$out". -fastqout "$out".Filtered.fastq
done

# "samples should be pooled before dereplication"
cat *Filtered.fastq >> allsamples.fastq

#dereplicate
usearch -fastx_uniques allsamples.fastq -fastqout allsamples.deRep.fastq -fastaout allsamples.deRep.fasta -sizeout -relabel zotu -minuniquesize 10

#run the unoise3 alogrithm
usearch -unoise3 allsamples.deRep.fasta -zotus zotus.fa -tabbedout unoise3.summary.txt

#run the otu tab command to get an otu table

usearch -otutab allsamples.fastq -zotus zotus.fa -otutabout Paua-eDNA.unoise.zotu.txt -mapout zmap.txt -notmatched Paua-eDNA.zotutab.unmapped.fasta


# a more fleshed out script
# Filters at 0.05 and then trims by 100bp / and then denoises

