#!/bin/bash -e
#SBATCH -A uoo02328
#SBATCH -J trim-emtDNA
#SBATCH --time 1:00:00
#SBATCH -N 1
#SBATCH -c 2
#SBATCH -n 1
#SBATCH --mem=16
#SBATCH --partition=large
#SBATCH --mail-user=clare.adams@postgrad.otago.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --output trim.%j.out # CHANGE map1 part each run
#SBATCH --error trim.%j.err # CHANGE map1 part each run


named="/nesi/nobackup/uoo02328/clare/Paua-e-mtDNA/CutAdapt4/AllSamples/"

cd "${named}" || exit

module load BBMap/38.81-gimkl-2020a

#mkdir BBmerged

#gunzip "${named}"*.gz

for sample in "${named}"*.1.fastq
do
        out=$(basename "$sample" .1.fastq)
        echo "$out"
        reformat.sh in1="${out}".1.fastq in2="${out}".2.fastq out1="${out}".filt1.fastq out2="${out}".filt2.fastq
        bbmerge.sh in1="${out}".filt1.fastq in2="${out}".filt2.fastq out=BBmerged/"${out}".merged.fastq outinsert=BBmerged/"${out}".stats.txt
done





module load cutadapt/2.10-gimkl-2020a-Python-3.8.2

cd /nesi/nobackup/uoo02328/clare/Paua-e-mtDNA/CutAdapt4/AllSamples/BBmerged || exit

for sample in *.merged.fastq
do
        out=$(basename "$sample" .merged.fastq)
        echo "$out"
        cutadapt -g CCTTGGTAACTAAACACACC -a CCTCCCTATTAGTCCTAAATC -m 268 -n 2 --cores=0 -o "${out}".trim.fastq "${out}".merged.fastq
done







cd /nesi/nobackup/uoo02328/clare/Paua-e-mtDNA/CutAdapt4/AllSamples/BBmerged/trimmed/2021-ee0.10 || exit

module load USEARCH/11.0.667-i86linux32

for sample in *.trim.fastq
do
        out=$(basename "$sample" .trim.fastq)
        echo "$out"
        usearch -fastq_filter "$sample" -fastq_maxee 0.1 -relabel "$out". -fastqout "$out".Filtered10.fastq
done

cat *Filtered10.fastq >> allsamples10.fastq

# now they will be dereplicated with a number
usearch -fastx_uniques allsamples10.fastq -fastqout allsamples10.deRep.fastq -fastaout allsamples10.deRep.fasta -sizeout -relabel Uniq -minuniquesize 10

# now they will be denoised
usearch -unoise3 allsamples10.deRep.fasta -zotus zotus10.fa -tabbedout unoise3.summary1.txt

#and a ZOTU table (ignore 97% clustering) will be made

#usearch -otutab allsamples.fastq -otus allsamples.deRep.fasta -otutabout NZFS-eDNA.otutab.txt -mapout map.txt -notmatched NZFS-eDNA.otutab.unmapped.fasta
usearch -otutab allsamples1.fastq -zotus zotus1.fa -otutabout Paua.eDNA.unoise.1.txt -mapout zmap.txt -notmatched Paua.eDNA.zotutab.unmapped1.fasta


#OTU ID Nothing5        Ohau1   Nothing6        Ohau2   Ohau3   Ohau5   Ohau7   Ohau8   Ohau9   Sample_10_1     Sample_10_2     Sample_10_3     Sample_1_1      Sample_1_2      
#Zotu1   1       3861    1       6818    7155    11444   1907    8143    2       0       6639    21      4484    33925   10820   8363    15614   12582   27649   19328   12
#Zotu2   0       2972    1       4167    10820   3057    10351   16677   255     232     8       4283    27535   1193    3468    4702    22645   14527   82      17193   131
#Zotu3   0       573     0       2641    4023    3571    12647   4948    52      0       1       1538    0       3       21628   0       5       309     0       16394   21
#Zotu4   0       1       0       0       224     7026    0       0       0       0       0       0       0       0       0       0       737     0       0       0       0
