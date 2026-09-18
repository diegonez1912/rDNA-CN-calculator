#!/usr/bin/env bash
#SBATCH --account=CancerEvolution
#SBATCH --time=12:00:00
#SBATCH --mem=64g
#SBATCH --cpus-per-task=16

# Load conda environment
set -euo pipefail
source /home/diegoonez/miniforge3/etc/profile.d/conda.sh
conda activate rdna_env

host_dir=/faststorage/project/CancerEvolution_shared/Projects/Diego
mkdir -p ${host_dir}/ref # Create a directory to store reference files (e.g. FASTA, GFF, BED) that will be used for mapping and analysis


#########################################################################
# Get the 45S rDNA reference
#########################################################################

# Download the reference file --> Canonical Human 45S rDNA repeat unit
  # This file (U13369.1) contains the full 45S rDNA repeat unit, which is approximately 43 kb in length.
  # The 45S rDNA repeat unit includes the 18S, 5.8S, and 28S rRNA genes, as well as the internal and external transcribed spacers.
echo "Downloading U13369.1 (45S rRNA) from GenBank..."
efetch -db nucleotide -id U13369.1 -format fasta > ${host_dir}/ref/U13369.1.fasta

# Index file, extract and concatenate regions
  # Due to the large size of the original 45S rDNA repeat unit (43 kb), we will create a modified reference that includes the first 14 kb and the last 2 kb of the sequence.
  # Thus, we will have a ~16 kb reference that captures the key regions of interest while being more manageable for downstream analyses.
echo "Concatenating to make rearranged reference..."
samtools faidx ${host_dir}/ref/U13369.1.fasta # This creates the .fai index file needed for region extraction
samtools faidx ${host_dir}/ref/U13369.1.fasta U13369.1:41021-42999 > ${host_dir}/ref/part1.fa # Extract the last 20 kb (41021-42999) of the original sequence
samtools faidx ${host_dir}/ref/U13369.1.fasta U13369.1:1-14000   > ${host_dir}/ref/part2.fa # Extract the first 14 kb (1-14000) of the original sequence

# Merge the sequences
{
  echo ">45S_rDNA_reference (U13369.1_41021-42999_1-14000)"
  cat "${host_dir}/ref/part1.fa" "${host_dir}/ref/part2.fa" | grep -v "^>" \
    | awk '
        { seq = seq $0 }
        END {
          for (i = 1; i <= length(seq); i += 70)
            print substr(seq, i, 70)
        }
      '
} > "${host_dir}/ref/45S_U13369.1_modified_16kb.fasta"

# Remove intermediate files
rm -f \
  "${host_dir}/ref/U13369.1.fasta" \
  "${host_dir}/ref/U13369.1.fasta.fai" \
  "${host_dir}/ref/part1.fa" \
  "${host_dir}/ref/part2.fa"

echo "45S rDNA reference built successfully:"
ls -lh ${host_dir}/ref/45S_U13369.1_modified_16kb.fasta

#########################################################################
# Get the 5S rDNA reference
#########################################################################

# Downloading reference file --> Canonical Human 5S rDNA repeat unit
  # This file (X12811.1) contains the canonical human 5S rDNA repeat unit, which is approximately 2 kb in length.
  # The 5S rDNA repeat unit includes the 5S rRNA gene and its associated regulatory regions.
echo "Downloading X12811.1 (5S rRNA) from GenBank..."
efetch -db nucleotide -id X12811.1 -format fasta > ${host_dir}/ref/5S_X12811.1.fasta

echo "5S rDNA reference downloaded successfully:"
ls -lh ${host_dir}/ref/5S_X12811.1.fasta

# Combine the two fasta files and index
  # We will combine the 5S and 45S rDNA reference sequences into a single FASTA file for easier mapping and analysis.
  # The BWA index will create several index files that are used for mapping, like: .rnd, .bwt, .pac, .ann, .amb
cat ${host_dir}/ref/5S_X12811.1.fasta ${host_dir}/ref/45S_U13369.1_modified_16kb.fasta > ${host_dir}/ref/rDNA_combined.fasta
bwa index ${host_dir}/ref/rDNA_combined.fasta


#########################################################################
# Download relevant reference files
#########################################################################

# Download the Ensembl reference genome (hg38)
  # This is the FASTA file for the human reference genome (GRCh38).
  # It will be used for mapping and analysis of the sequencing data.
wget -O ${host_dir}/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
  https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip ${host_dir}/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
samtools faidx ${host_dir}/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa # This creates the .fai index file needed for region extraction

# Download the genome annotation (GFF)
  # This is the GTF file for the human reference genome (GRCh38).
  # It contains gene annotations and is used for telling us what biological features are in each genomic region.
wget -O ${host_dir}/ref/GCF_000001405.40_GRCh38.p14_genomic.gff.gz \
  https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gff.gz
gunzip ${host_dir}/ref/GCF_000001405.40_GRCh38.p14_genomic.gff.gz

# Download repetitive sequences from RepeatMasker (converted to BED file)
  # This is a checking point for DNA quality. 
  # It screens DNA sequences for interspersed repeats and low complexity DNA sequences.
  # In other words, it tells us: “Which parts of chromosome 1 are genes, transcripts, exons, CDS, etc.?”
wget -O ${host_dir}/ref/rmsk_hg38.txt.gz \
  http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz
zcat ${host_dir}/ref/rmsk_hg38.txt.gz | awk -v OFS='\t' '{sub(/^chr/, "", $6); print $6, $7, $8, $12}' > ${host_dir}/ref/repeatmasker_hg38.bed

rm ${host_dir}/ref/rmsk_hg38.txt.gz


#########################################################################
# Get exonic and intronic regions
#########################################################################

python scripts/build_background_ref.py \
  --genome ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --gff ref/GCF_000001405.40_GRCh38.p14_genomic.gff \
  --rmsk ref/repeatmasker_hg38.bed \
  --outdir /faststorage/project/CancerEvolution_shared/Projects/Diego/ref

cat ${host_dir}/ref/bg_chr1_exon.bed ${host_dir}/ref/bg_chr1_intron.bed > ${host_dir}/ref/bg_chr1_exons_introns.bed
cat ${host_dir}/ref/bg_acro_exon.bed ${host_dir}/ref/bg_acro_intron.bed > ${host_dir}/ref/bg_acro_exons_introns.bed
