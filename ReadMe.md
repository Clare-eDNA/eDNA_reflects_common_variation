# eDNA reflects common variation

Please cite the paper if any of this was useful to you.

Please consider using the files in the order of: 

CutAdapt_F/R_Barcodes - to help demultiplex the data
2-4_data_analyzing_final_script.sh - to demultiplex, merge, filter, and put the data into an OTU table
Note that 1_CutAdaptDemux.sh and 3_filt+unoise.sh scripts can be modified and combined for the above script

DataMunge-eDNAamplicon.R reads in the OTU table (Paua.eDNA.ee10.csv) and makes some changes for readability and labeling.
eDNA-Paua-Graphs.R helps to create graphs for the paper.

Haplotype Accumulation graph and Maps are stand-alone.
