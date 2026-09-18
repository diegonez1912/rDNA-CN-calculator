#!/usr/bin/env bash
#SBATCH --account=
#SBATCH --time=50:00:00
#SBATCH --mem=32g
#SBATCH --cpus-per-task=26


# Load conda environment

THREADS=8 # $SLURM_CPUS_PER_TASK

########################################################################
# DEFINE INPUT AND OUPUT ARGUMENTS, AS WELL AS PATHS TO REFERENCE FILES
########################################################################

# Define input arguments
BAM="$1" # First command line argument (BAM file path of the sample to process)
PROJECT_ID="$2" # Second command line argument (project name, e.g. GBM_GL) used to organize output directories

# Define output configuration
REF_DIR="ref" # Directory where reference files (FASTA, BED) are stored
# LOGS="logs"
OUT_DIR="output/${PROJECT_ID}" # Directory to store rDNA depth outputs (e.g. position-wise depth files)

# Create directories (if they are needed)
# mkdir -p "${LOGS}"
mkdir -p "${OUT_DIR}"

# Sample name extraction from BAM filename
BAM_BASENAME="$(basename "${BAM}")" # Extract the base name of the BAM file (e.g. TCGA-4Z-AA7S-10A.bam)
SAMPLE="${BAM_BASENAME%%.*}" # Remove the extension to get the sample name (e.g. TCGA-4Z-AA7S-10A)
# OUT_DIR="${OUT_DIR}/${SAMPLE}"
# mkdir -p "${OUT_DIR}"

# # Redirect stdout and stderr to log
# LOG="${LOGS}/${SAMPLE}.log"
# exec > >(tee -a "$LOG") 2>&1

# Reference files
RDNA_REF="${REF_DIR}/rDNA_combined.fasta" # Combined FASTA with 5S and 45S rDNA reference sequences (from get_reference_seq.sh)
BG_CHR1_BED="ref/bg_chr1_exons_introns.bed" # Single-copy regions on chr1 used for 5S CN estimation (from build_background_ref.py)
BG_OTHER_BED="ref/bg_acro_exons_introns.bed" # Single-copy regions on chr13/14/15/21/22 used for 45S CN estimation (from build_background_ref.py)


########################################################################
# SLICE rDNA REGIONS AND MAP TO rDNA REFERENCE
########################################################################

# rDNA slice coordinates in the *genome-aligned* BAM (as in the paper)
RDNASCAF_45S="chrUn_GL000220v1" # 45S reads mostly here
FIVE_S_REGION="chr1:226743523-231781906" # 1q42 +/- 2 Mb around 5S locus

# Names of rDNA references inside rDNA_combined.fasta
RDNA_5S_CONTIG="X12811.1"
RDNA_45S_CONTIG="45S_rDNA_reference"

# Check if output already exists
FINAL_5S="${OUT_DIR}/${SAMPLE}.5S.BRD_norm_depth.tsv"
FINAL_45S="${OUT_DIR}/${SAMPLE}.45S.BRD_norm_depth.tsv"

if [[ -f "${FINAL_5S}" && -f "${FINAL_45S}" ]]; then
  echo "[INFO] Final output already exists for ${SAMPLE}. Skipping."
  exit 0
fi

# Check if all reference files and their paths are there
for f in "${BAM}" "${RDNA_REF}" "${BG_CHR1_BED}" "${BG_OTHER_BED}"; do
  if [[ ! -f "$f" ]]; then
    echo "ERROR: Required file not found: $f" >&2
    exit 1
  fi
done

# Check if the bam file is indexed --> If not found, create the index file (.bai)
if [[ ! -f "${BAM%.bam}.bai" && ! -f "${BAM}.bai" ]]; then 
  samtools index "${BAM}"; 
fi

echo "========================================"
echo "[START] $(date)"
echo "[INFO] Processing sample: ${SAMPLE}"
echo "[INFO] Input BAM: ${BAM}"
echo "========================================"

###############################
# Step 1: Slice rDNA regions from the genome-aligned BAM file (the sample) and convert to FASTQ
###############################

# Output FASTQ file for the sliced reads that will be mapped to the rDNA reference
SAMPLE_FQ="${OUT_DIR}/${SAMPLE}.rDNA_slice.fastq"
echo "[INFO] Slicing BAM"

# Extract alignments from the BAM file that overlap with the 45S rDNA scaffold and the 5S region on chr1
  # --fetch-pairs --> If one read of a pair overlaps with the specified scaffold/region, both reads of the pair will be included in the output
  # samtools fasq --> Convert extracted alignments to FASTQ format
  # SAMPLE_FQ --> File containing FASTQ reads OF THE WHOLE GENOME that aligned to the 45S rDNA scaffold or the 5S region on chr1
    # They are just sequences that might appear in different parts of the genome (different copies of rDNA), so we need to map them to one rDNA reference (step 2)
samtools view --fetch-pairs  -b "${BAM}" "$RDNASCAF_45S" "$FIVE_S_REGION" \
  | samtools fastq -N -o "${SAMPLE_FQ}"


########################################
# Step 2: Map rDNA-sliced reads to rDNA reference to get position-wise depth
########################################

# Output BAM file for the alignments of the sliced rDNA reads to the rDNA reference
  # samtools depth only works on aligned BAM files
RDNA_BAM="${OUT_DIR}/${SAMPLE}.rDNA_onRef.sorted.bam"
echo "[INFO] Mapping FASTQ to rDNA reference"

# Map sliced FASTQ rDNA reads to a single-copy rDNA reference, so that we see how many rDNA reads support that canonical rDNA reference sequence
  # By doing this, the coverage of rDNA becomes proportional to the copy number of rDNA in the sample (even if they come from different copies of the same rDNA reference)
  # bwa mem --> Map the FASTQ reads to the rDNA reference using BWA-MEM algorithm
  # -p --> Indicates that the input is interleaved paired-end FASTQ (both mates in one file)
  # samtools sort --> Sort the resulting alignments by coordinate and output as BAM
  # RDNA_BAM --> BAM file containing the alignments of the sliced rDNA reads to the rDNA reference, sorted by coordinates
bwa mem -t "$THREADS" -p ${RDNA_REF} ${SAMPLE_FQ} \
  | samtools sort -o "${RDNA_BAM}"
samtools index "${RDNA_BAM}"


########################################
# Step 3. Estimate BRD for single-copy regions
########################################

echo "[INFO] Estimating background read depth (BRD) for single-copy regions"

# BRD_CHR1 --> Average read depth in single-copy regions on chr1 (used for 5S CN estimation)
  # samtools depth -b --> Compute read depth at positions specified in the BED file (BG_CHR1_BED) for the original BAM file (the sample)
  # awk --> Calculate the average BRD across all positions/base-pairs by summing the depth values and dividing by the number of base-pairs in the reference-region
BRD_CHR1=$(samtools depth -b "${BG_CHR1_BED}" "${BAM}" \
  | awk '{sum+=$3; n++} END{print sum/n}')

# BRD_OTHER --> Average read depth in single-copy regions on chr13/14/15/21/22 (used for 45S CN estimation)
BRD_OTHER=$(samtools depth -b "${BG_OTHER_BED}" "${BAM}" \
  | awk '{sum+=$3; n++} END{print sum/n}')

echo "[INFO] BRD chr1 (5S) : ${BRD_CHR1}"
echo "[INFO] BRD chr13+ (45S): ${BRD_OTHER}"


########################################
# Step 4: Get position wise BRD for rDNA aligned reads (from step 2)
########################################

echo "[INFO] Computing normalized position-wise depth on rDNA samples"

# 5S --> Normalize by BRD from chr1 single-copy regions
  # samtools depth -a -r --> Compute read depth at all positions specified by the region (RDNA_5S_CONTIG) for the BAM file (RDNA_BAM)
    # As samtools depth -a includes zeros, this mean includes uncovered positions
  # awk --> Normalize depth at each position by dividing the depth value with the mean BRD value from single-copy regions on chr1 (BRD_CHR1)
samtools depth -a -r "${RDNA_5S_CONTIG}" "${RDNA_BAM}" \
  | awk -v brd="${BRD_CHR1}" 'BEGIN{OFS="\t"} {print $1,$2,$3/brd}' \
  > "${OUT_DIR}/${SAMPLE}.5S.BRD_norm_depth.tsv"

# 45S --> Normalize by BRD from chr13/14/15/21/22 single-copy regions
samtools depth -a -r "${RDNA_45S_CONTIG}" "${RDNA_BAM}" \
  | awk -v brd="${BRD_OTHER}" 'BEGIN{OFS="\t"} {print $1,$2,$3/brd}' \
  > "${OUT_DIR}/${SAMPLE}.45S.BRD_norm_depth.tsv"

# Cleanup intermediate files
rm -f "${SAMPLE_FQ}"
rm -f "${RDNA_BAM}" "${RDNA_BAM}.bai"

echo "[DONE]"



