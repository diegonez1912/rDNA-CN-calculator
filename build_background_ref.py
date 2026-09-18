#!/usr/bin/env python3
import argparse, os, subprocess
import pyranges as pr
import pandas as pd
import pysam
import numpy as np

"""
Key libraries: 
- argparse for command-line argument parsing
- os and subprocess for file handling and running external commands (makeblastdb, blastn)
- pyranges for genomic intervals manipulation
- pandas for data handling
- pysam for FASTA handling (fetch sequences from FASTA by coordinates)
- numpy for efficient numeric operations
"""

"""
Run for hg38 as:
scripts/build_background_ref.py 
--genome ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa 
--gff ref/GCF_000001405.40_GRCh38.p14_genomic.gff 
--rmsk ref/repeatmasker_hg38.bed 
--outdir /faststorage/project/CancerEvolution_shared/Projects/Diego/scripts/ref
"""

##############################################################
# Functions
##############################################################


def compute_introns_from_transcripts_minus_exons(gff_gr, exons_all):
    """
    GOAL --> Construct a set of intronic intervals for each gene
        - Introns(transcript) = Transcript_span - exons(transcript)
    
    INPUTS:
        - gff_gr --> PyRanges dataframe with at least rows for mRNA and exon features, with columns including Chromosome, Start, End, Strand, Feature, gene, ID, Parent
        - exons_all --> PyRanges dataframe containing exons (Feature == "exon") // Includes a "gene" column and Parent (transcript ID)
    
    OUTPUT --> PyRanges dataframe with intron intervals, with an added "gene" column
    """

    # Filter transcripts --> Keep only mRNA molecules (the ones with coding regions) // .copy() to avoid mutating the original gff_gr
    tx = gff_gr[gff_gr.Feature == "mRNA"].copy()

    results = [] # for collecting introns per gene
    cols = ["Chromosome","Start","End","Strand"]

    # Loop through genes in the df and their ID transcripts:
        # Convert it to pandas (.df) for easier manipulation and group mRNA transcripts by gene
        # "gene_sub" --> All mRNA rows for that specific gene (introns and exons)
    for gene, gene_sub in tx.df.groupby("gene"):
        
        # Pull all exon rows for that specific gene from "exons_all" --> "ex_gene_sub" (grouped by gene)
        ex_gene_sub = pr.PyRanges(exons_all.df[exons_all.df["gene"] == gene])
        if ex_gene_sub.empty:
            continue

        gene_introns = []

        # Loop through transcript IDs of that specific gene:
            # Exons have "Parent" equal to that transcript ID
            # trans = a transcript ID (e.g. NM_... style)
            # trans_sub = All mRNA row(s) in gene_sub for that transcript ID (usually is only one row/mRNA per transcript ID)
        for trans, trans_sub in gene_sub.groupby('ID'):
            # Get exons of the transcript ID (we might not use all exons of the gene due to alternative splicing)
            sub_ex = ex_gene_sub[ex_gene_sub.Parent == trans]
            if sub_ex.empty:
                continue
            
            # Compute introns --> A single gene can have multiple transcripts with its own set of exons:
                # That's why we compute introns per ID transcript and not only per gene
                # .substract() --> PyRanges method to subtract intervals (returns parts of the script not covered by exons)
            intron = pr.PyRanges(trans_sub[cols]).subtract(sub_ex[cols])

            # Get exons from other transcripts, and remove introns that overlap other exons
                # overlap(..., invert = True) --> Keep only regions that do NOT overlap the provided intervals
            other_trans = ex_gene_sub[ex_gene_sub.Parent != trans]
            if not other_trans.empty:
                intron = intron.overlap(other_trans[cols], invert = True)
            # Add the intron to the list
            if not intron.empty:
                gene_introns.append(intron)

        if gene_introns:
            gene_introns = pr.concat(gene_introns).merge() # pr.concat --> Concatenates all intron PyRanges from isoforms
            gene_introns = gene_introns.df # convert to dataframe to add "gene" column
            gene_introns["gene"] = gene # add "gene" column
            results.append(gene_introns)

    # Concatenate all genes’ introns into one DataFrame, wrap in PyRanges, and sort
    return pr.PyRanges(pd.concat(results, ignore_index = True)).sort()


def pick_largest_exon_per_overlap(gr_exon):
    """
    GOAL --> For each gene, pick the largest exon among overlapping exons from different transcripts (isoforms)
    
    INPUT --> PyRanges dataframe with feature == 'exon' // has "gene" column

    OUTPUT --> PyRanges of selected exons (one per overlapping cluster per gene)
    """
    
    # Cluster overlapping exons per gene // We also keep gene_id in the clustering key.
    clusters = gr_exon.cluster(by = "gene")  # adds a column 'Cluster' that assigns a cluster ID to overlapping exons within the same gene
    df = clusters.df.copy() # Convert "clusters" to dataframe for easier manipulation

    # Keep the largest exon within each (gene, Cluster)
    df["len"] = (df["End"] - df["Start"]).astype(int) # Compute exon length
    idx = df.groupby(["gene", "Cluster"])["len"].idxmax() # Find the index (exon) with the maximum length for each (gene, Cluster)
    picked = pr.PyRanges(df.loc[idx, ["Chromosome", "Start", "End", "Strand", "gene"]].sort_values(["Chromosome", "Start", "End"]))

    return picked


def load_repeatmasker_bed(rmsk_path):
    """
    GOAL --> Load RepeatMasker BED file and return a PyRanges object with merged intervals (overlapping repeats merged into one interval)
        - We use this to ensure that the candidate background regions behave like single-copy DNA regions
        - RepeatMasker is a annotation file of repetitive elements in the genome (transposable elements, simple repeat, low-complexity regions, etc.)

    INPUT --> Path to RepeatMasker BED file (with at least columns: chrom, start, end)

    OUTPUT --> PyRanges object with merged intervals of repeats
    """
    
    df = pd.read_csv(rmsk_path, 
                     sep="\t", 
                     header=None,
                     usecols=[0, 1, 2], # keep only 3 first columns (chrom, start, end)
                     dtype={0: "string", 1: np.int32, 2: np.int32}, # specify dtypes for memory efficiency
                     names=["Chromosome", "Start", "End"])
    
    return pr.PyRanges(df).merge()


def load_gff_and_get_coding_transcripts(gff_path):
    """
    GOAL --> Load the GFF annotation file and produce a Pyranges file containing protein-coding genes, their mRNA transcripts and exon features

    INPUT --> Path to GFF file (with standard GFF columns and attributes including gene_biotype, ID, Parent, gene)

    OUTPUT --> PyRanges dataframe  with columns including Chromosome, Start, End, Strand, ID, Parent, gene
    """

    # Map RefSeq accessions to chromosome numbers (e.g. NC_000001.11 --> 1) for easier handling    
    refseq_to_chr = {
        "NC_000001.11": "1", "NC_000002.12": "2", "NC_000003.12": "3", "NC_000004.12": "4",
        "NC_000005.10": "5", "NC_000006.12": "6", "NC_000007.14": "7", "NC_000008.11": "8",
        "NC_000009.12": "9", "NC_000010.11": "10", "NC_000011.10": "11", "NC_000012.12": "12",
        "NC_000013.11": "13", "NC_000014.9": "14", "NC_000015.10": "15", "NC_000016.10": "16",
        "NC_000017.11": "17", "NC_000018.10": "18", "NC_000019.10": "19", "NC_000020.11": "20",
        "NC_000021.9": "21", "NC_000022.11": "22", "NC_000023.11": "X", "NC_000024.10": "Y"
    }

    # Load the GFF file as a PyRanges object, and convert to pandas DataFrame for easier manipulation
    gff_df = pr.read_gff3(gff_path, as_df = True)

    # Rename RefSeq accessions to chromosome numbers
    gff_df["Chromosome"] = gff_df["Chromosome"].replace(refseq_to_chr)

    # Get a list of protein coding genes --> Filter to rows annotated as protein_coding genes
    genes_pc = gff_df[gff_df.gene_biotype == "protein_coding"]
    gene_ids = set(genes_pc["ID"].dropna()) # get the set of gene IDs for protein coding genes (dropna to avoid NaN values)

    # Get all rows annotated as mRNA transcripts and keep only those whose parent gene is a protein coding gene (gene_id = "Parent") 
    tx = gff_df[gff_df.Feature == "mRNA"].copy()
    tx = tx[tx["Parent"].isin(gene_ids)]
    tx_ids = set(tx["ID"].dropna()) # Collect transcript IDs of mRNA transcripts

    # Make a final df with mRNA and exons:
        # Start from all mRNA and exon rows
        # Filter to mRNA or exon of the selected mRNAs that are children of the selected protein coding genes
    tx_exons = gff_df[gff_df.Feature.isin(["mRNA", "exon"])].copy() 
    tx_exons = tx_exons[(tx_exons['ID'].isin(tx_ids)) | (tx_exons['Parent'].isin(tx_ids))]

    # Reduce to relevant columns and sort
    tx_exons = tx_exons[['Chromosome', 'Feature', 'Start', 'End', 'Strand', 'ID', 'Parent', 'gene']].copy()

    # Convert to PyRanges for interval operations and return
    return pr.PyRanges(tx_exons)


def build_candidates(gff_path, rmsk_path):
    """
    GOAL --> Build 4 sets of candidate regions for background reference:
    1) Exons on chromosome 1
    2) Introns on chromosome 1
    3) Exons on acrocentric chromosomes (13, 14, 15, 21, 22)
    4) Introns on acrocentric chromosomes (13, 14, 15, 21, 22)

    INPUTS:
        - gff_path --> Path to GFF annotation file for GRCh38 (with attributes including gene_biotype, ID, Parent, gene)
        - rmsk_path --> Path to RepeatMasker BED file for hg38 (with at least columns: chrom, start, end)

    OUTPUTS:
        - cands_chr1_exon --> PyRanges of candidate exons on chromosome 1
        - cands_chr1_intron --> PyRanges of candidate introns on chromosome 1
        - cands_acro_exon --> PyRanges of candidate exons on acrocentric chromosomes
        - cands_acro_intron --> PyRanges of candidate introns on acrocentric chromosomes
    """

    # Load the gff with genome annotation of Ensembl
    gff_gr = load_gff_and_get_coding_transcripts(gff_path)

    # Extract exons and introns on chromosomes 1, 13, 14, 15, 21 and 22
    gff_gr = gff_gr[gff_gr.Chromosome.isin(["1", "13","14","15","21","22"])]

    # Extract exons of protein coding transcripts (mRNAs)
    exons_all = gff_gr[gff_gr.Feature == "exon"]

    # Pick largest exon in overlapping isoforms (by gene + cluster)
    exons = pick_largest_exon_per_overlap(exons_all)

    # Extract introns from transcript spans minus union of exons
        # Remove intron pieces that are exonic in other isoforms of the same gene
    introns = compute_introns_from_transcripts_minus_exons(gff_gr, exons_all)

    # Split exons and introns of reference candidates in chromosome 1
    ex_chr1 = exons[exons.Chromosome == '1']
    in_chr1 = introns[introns.Chromosome == '1']

    # Split exons and introns of reference candidates in acrocentric chromosomes (13, 14, 15, 21, 22)
    acro_chroms = ["13","14","15","21","22"]
    ex_acro = exons[exons.Chromosome.isin(acro_chroms)]
    in_acro = introns[introns.Chromosome.isin(acro_chroms)]

    # Load repeats, so that we ensure that reference candidates are single-copy DNA regions
    rmsk_gr = load_repeatmasker_bed(rmsk_path)

    # Substract repeats from reference candidates:
        # overlap(..., invert = True) --> Removes anything that overlaps repeats
        # .merge() --> Unions overlapping/adjacent resulting intervals
    cands_chr1_exon, cands_chr1_intron = ex_chr1.overlap(rmsk_gr, invert = True).merge(), in_chr1.overlap(rmsk_gr, invert = True).merge()
    cands_acro_exon, cands_acro_intron = ex_acro.overlap(rmsk_gr, invert = True).merge(), in_acro.overlap(rmsk_gr, invert = True).merge()

    return cands_chr1_exon, cands_chr1_intron, cands_acro_exon, cands_acro_intron


def len_filter(gr, minlen = 300, maxlen = 10000):
    """
    GOAL --> Filter PyRanges to keep only reference candidates with length between minlen and maxlen (inclusive)
    """
    
    df = gr.df.copy() # Convert to dataframe for easier manipulation
    df["len"] = (df["End"] - df["Start"]).astype(int) # Compute interval length
    df = df[(df["len"] >= minlen) & (df["len"] <= maxlen)] # Filter to intervals with length between minlen and maxlen
    
    return pr.PyRanges(df[["Chromosome","Start","End"]]).sort() # Convert back to PyRanges and sort


def trim(gr, trim_bp = 50):
    """
    GOAL --> Trim a specified number of base pairs (trim_bp) from both ends of each interval in the PyRanges
        - Filter out intervals that become invalid (end <= start) after trimming
        - Boundary regions are more likely to have alignment issues, so we trim them to ensure more reliable single-copy regions
    """
    
    df = gr.df.copy() # Convert to dataframe for easier manipulation
    df["Start"] = df["Start"].astype(int) + trim_bp # Trim from start
    df["End"] = df["End"].astype(int) - trim_bp # Trim from end
    df = df[df["End"] > df["Start"]] # Filter out intervals that become invalid after trimming (end must be greater than start)
    
    return pr.PyRanges(df[["Chromosome","Start","End"]]).sort() # Convert back to PyRanges and sort


def write_fasta_from_pyranges(gr, genome_fa, out_fa):
    """
    GOAL --> Extract DNA sequences for each reference candidate in build_candidates(), by using the reference genome (FASTA)
        - After extracting the sequences, write them to a FASTA file where each record header is in the format 'chrom:start-end' 
    
    INPUTS:
        - gr --> PyRanges dataframe of reference candidates // Columns: Chromosome, Start, End (1-based coordinates)
        - genome_fa --> Path to the reference genome FASTA file (indexed with .fai)
        - out_fa --> Path to the output FASTA file to write the sequences of the reference candidates
    
    OUTPUT:
        - A FASTA file at out_fa containing the sequences of the reference candidates, with headers in the format 'chrom:start-end'
    """
    
    fa = pysam.FastaFile(genome_fa) # Open reference genome FASTA file with pysam (requires .fai index)
    
    with open(out_fa, "w") as out: # Open output FASTA file for writing
        for _, r in gr.df.iterrows(): # Iterate through each row of the PyRanges dataframe (reference candidates)
            
            # Extract chromosome, start, and end from the row (convert chromosome to string for pysam)
            chrom = str(r.Chromosome)
            start = int(r.Start)
            end   = int(r.End)
            
            # Fetch the sequence from the reference genome using 0-based coordinates (pysam uses 0-based coordinates, while PyRanges uses 1-based)
            seq = fa.fetch(chrom, start, end)
            # Write the sequence to the output FASTA file with a header in the format 'chrom:start-end'
            out.write(f">{chrom}:{start}-{end}\n")
            
            # Wrap sequence at 60 base-pairs per line
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + "\n")
    
    # Close the FASTA file
    fa.close()


def run_makeblastdb(genome_fa, db_prefix):
    """
    GOAL --> Create a BLAST nucleotide database from the reference candidate FASTA file
        - By creating a BLAST database, we can blast candidate sequences against the whole genome to check if it has multiple alignments (i.e. if it's not single-copy)

    INPUTS:
        - genome_fa --> Path to the reference candidate FASTA file (indexed with .fai)
        - db_prefix --> Prefix for the output BLAST database files (e.g. "hg38_blastdb/hg38")

    OUTPUT --> BLAST database files with the specified prefix (e.g. "hg38_blastdb/hg38.nhr", "hg38_blastdb/hg38.nin", "hg38_blastdb/hg38.nsq")
    """
    
    # Ensure the output directory for the BLAST database exists
    os.makedirs(os.path.dirname(db_prefix), exist_ok = True)
    
    subprocess.run([
        "makeblastdb", # Command-line tool to create a BLAST database
        "-in", genome_fa, # Input FASTA file (reference genome)
        "-dbtype", "nucl", # Specify that the database is for nucleotide sequences
        "-parse_seqids", # Parse sequence IDs from the FASTA headers (required for blastn to work with sequence IDs)
        "-out", db_prefix # Output prefix for the BLAST database files
    ], check = True)


def run_blastn(query_fa, db_prefix, out_tsv, threads = 8, evalue = 1e-6):
    """
    GOAL --> Run BLASTN to align the candidate sequences (query_fa) against the reference genome BLAST database (db_prefix) and calculate BLAST hits

    INPUTS:
        - query_fa --> Path to the FASTA file containing the candidate sequences to be aligned (from write_fasta_from_pyranges())
        - db_prefix --> Prefix of the BLAST database to align against (from run_makeblastdb())
        - out_tsv --> Path to the output TSV file where the BLAST results will be written
        - threads --> Number of threads to use for BLAST (default: 8)
        - evalue --> E-value threshold for reporting BLAST hits (default: 1e-6)

    OUTPUT --> TSV file with BLAST results (columns: qseqid sseqid pident length evalue bitscore qstart qend sstart send)
        - qseqid --> Query ID (your chrom:start-end)
        - sseqid --> Subject ID (genome contig)
        - Identity, alignment length, evalue, bitscore, and coordinates
    """
    
    subprocess.run([
        "blastn",
        "-query", query_fa,
        "-db", db_prefix,
        "-out", out_tsv,
        "-outfmt", "6 qseqid sseqid pident length evalue bitscore qstart qend sstart send", # Output format: tab-delimited with specified columns
        "-evalue", str(evalue),
        "-num_threads", str(threads),
        "-task", "megablast"
    ], check = True)


def filter_single_copy_from_blast(gr, blast_tsv):
    """
    GOAL --> Filter the reference candidates to keep only those that have a single BLAST hit
        - One hit --> Assume that they are single-copy regions (align only to themselves in the genome)
        - Multiple hits --> Likely not single-copy, so we discard them
    
    INPUTS:
        - gr --> PyRanges dataframe of reference candidates (before filtering)
        - blast_tsv --> Path to the BLAST output TSV file (from run_blastn())

    OUTPUT --> PyRanges dataframe of reference candidates that have only one BLAST hit
    """
    
    # Load the BLAST output file into a dataframe with appropriate column names
    df = pd.read_csv(blast_tsv, sep = "\t", header = None,
                     names = ["qseqid","sseqid","pident","length","evalue","bitscore","qstart","qend","sstart","send"])
    
    # Count the number of BLAST hits for each query sequence (qseqid), and keep only those with a single hit
    counts = df.groupby("qseqid").size()
    keep_ids = set(counts[counts == 1].index)

    # Reconstruct the original Pyranges dataframe with a new "id" column that matches the qseqid format (chrom:start-end)
        # Filter the original Pyranges to keep only rows whose "id" is in the set of single-hit qseqids
    gdf = gr.df.copy()
    gdf["id"] = gdf["Chromosome"].astype(str) + ":" + gdf["Start"].astype(str) + "-" + gdf["End"].astype(str)
    gdf = gdf[gdf["id"].isin(keep_ids)]

    # Return the filtered candidates as a PyRanges dataframe, sorted by chromosome and start position
    return pr.PyRanges(gdf[["Chromosome","Start","End"]]).sort()


##############################################################
# Code
##############################################################

def main():
    """
    GOAL --> Build a set of single-copy background reference regions for hg38, by following these steps:
        1) Build exon/intron candidates and subtract repeats (from build_candidates())
        2) Run makeblastdb (once per genome path) to create a BLAST database from the reference genome FASTA (from run_makeblastdb())
        3) For each possible reference candidate set (chr1 exon, chr1 intron, acro exon, acro intron):
            - Write a FASTA file as input for BLAST (headers: chrom:start-end)
            - Run BLAST to align candidates against the whole genome
            - Filter candidates to keep only those with a single BLAST hit (i.e. likely single-copy regions)
            - Filter candidates by length (e.g. between 300 and 10000 bp) and trim boundary regions (e.g. 50 bp from each end)
            - Write final candidates to BED file (columns: chrom, start, end)
    """
        
    # Define command-line arguments for the script
    parser = argparse.ArgumentParser(description="Build hg38 background regions (single-copy exons/introns)")
    parser.add_argument("--genome", required = True, help = "GRCh38 FASTA (indexed; *.fai present)")
    parser.add_argument("--gff",    required = True, help = "Annotation GFF for GRCh38")
    parser.add_argument("--rmsk",   required = True, help = "RepeatMasker BED for hg38")
    parser.add_argument("--outdir", required = True, help = "Output directory")
    args = parser.parse_args()

    # Make sure the output directory exists
    os.makedirs(args.outdir, exist_ok = True)

    # 1) Build exon/intron candidates and subtract repeats
    cands_chr1_exon, cands_chr1_intron, cands_acro_exon, cands_acro_intron = build_candidates(args.gff, args.rmsk)

    # 2) Run makeblastdb (once per genome path) to create a BLAST database
        # Set database path prefix (e.g. "hg38_blastdb/hg38")
        # Check if the BLAST database files already exist to avoid redundant computation
    db_prefix = os.path.join(args.outdir, "hg38_blastdb", "hg38")
    if not (os.path.exists(db_prefix + ".nhr") or os.path.exists(db_prefix + ".nin") or os.path.exists(db_prefix + ".nsq")):
        run_makeblastdb(args.genome, db_prefix)

    # 3) Name each set of reference candidates
    candidates = {
        "chr1_exon":  cands_chr1_exon,
        "chr1_intron": cands_chr1_intron,
        "acro_exon":   cands_acro_exon,
        "acro_intron": cands_acro_intron,
    }
    
    # Loop through each set of candidates and: run BLAST, filter by single-copy, length, and trim boundaries, and write final BED files
    for name, gr in candidates.items():
        
        # Define intermediate and final output file paths for this candidate set
        out_fa = f"{args.outdir}/{name}_candidates.fa"
        out_tsv = f"{args.outdir}/{name}_candidates.tsv"
        out_bed = f"{args.outdir}/bg_{name}.bed"

        # Run blast
        write_fasta_from_pyranges(gr, args.genome, out_fa)
        run_blastn(out_fa, db_prefix, out_tsv)

        # Filter by single-copy, length, and trim boundaries
        bg_untrim = filter_single_copy_from_blast(gr, out_tsv)
        filter_df = len_filter(bg_untrim)
        filter_trim_df = trim(filter_df)

        # Write final candidates to BED file (columns: chrom, start, end)
        final_df = filter_trim_df.df[["Chromosome","Start","End"]]
        final_df["Chromosome"] = 'chr' + final_df["Chromosome"].astype(str)
        final_df.to_csv(out_bed, sep = "\t", header = False, index = False)

        # Clean up intermediate files
        os.remove(out_fa)
        os.remove(out_tsv)

        print(f"[ok] Wrote:\n  {out_bed}")


if __name__ == "__main__":
    main()


