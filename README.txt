This repository contains:

## 1 ##
get_reference_seq.sh (Selecting reference rDNA and background sequences):
  1. Downloads and modifies the 45S and 5S rDNA sequences.
  2. Downloads reference data.
  3. Calls the script build_background_ref.py to retreive the exonic and intronic regions to compare to.

Run as: sbatch get_reference_seq.sh
- Submit the file from the project folder
- Submit as --> sbatch get_reference_seq.sh
- Outputs --> ls -lh ref/rDNA_combined.fasta // ls -lh ref/rDNA_combined.fasta.*


## 2 ##
calc_rDNA_depth.sh (Obtaining rDNA array reads, and calculate position wise coverage of 5S and 45S):
  1. Slice rDNA regions from bam file and output to fastq
  2. Map rDNA-sliced reads to rDNA reference sequences
  3. Estimate background read depth (BRD) for single-copy regions
  4. Get position wise depth for the aligned files normalized by BRD

Run as: sbatch calc_rDNA_depth.sh bam project


## 3 ##
collect_rDNA_CN.py (Gets the rDNA copy number for the rDNA subunits for all patients in the folder)
Run as: python scripts/collect_rDNA_CN.py --folder folder/ --output output_file.tsv

To run the whole thing, run:
1. get_reference_seq.sh (only once)
2. calc_rDNA_depth.sh
3. collect_rDNA_CN.py

