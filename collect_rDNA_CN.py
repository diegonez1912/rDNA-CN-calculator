#!/usr/bin/env python3
import argparse, os, subprocess
from pathlib import Path
import pandas as pd
import numpy as np

"""
Key libraries:
- argparse: For parsing command-line arguments
- os, subprocess: For interacting with the operating system and running shell commands
- pathlib: For handling filesystem paths in an object-oriented way
- pandas: For data manipulation and analysis
- numpy: For numerical operations, especially with arrays
"""

"""
Run as:
python scripts/collect_rDNA_CN.py 
--folder output/TCGA_PRAD/ 
--output output/summary_files/TCGA_PRAD_rDNA_copy_number.tsv
"""

def annotate_45s_unit(pos: pd.Series) -> pd.Series:
    """
    GOAL --> Given a Series of positions along the 45S reference, label each position as 18S, 5.8S, or 28S
    Returns a categorical Series with values in {18S, 5.8S, 28S} or NaN.
    """
    LABELS = ["18S", "5.8S", "28S"]

    # Coordinates for different subunits --> They correspond to positions on the modified 16kb 45S reference
    conditions = [
        pos.between(5636, 7506), # 18S
        pos.between(8602, 8758), # 5.8S
        pos.between(9914, 14948), # 28S
    ]
    
    # Return as a pandas series with the same index as input
    return pd.Series(np.select(conditions, LABELS, default=None), index=pos.index)


def collect_rDNA_CN(project_folder):

    # Define patient IDs:
        # We assume that the project folder contains files named like {patient_id}.5S.BRD_norm_depth.tsv and {patient_id}.45S.BRD_norm_depth.tsv
        # We can extract patient IDs by listing all .tsv files and taking the part before 
        # Path(...).glob("*.tsv") --> Lists all TSV files in the folder // f.name --> Filename only
        # .split(".")[0] --> Takes the part before the first dot as sample ID // Example: {patient_id}.5S.BRD_norm_depth.tsv --> {patient_id}
    patient_ids = list({f.name.split(".")[0] for f in Path(project_folder).glob("*.tsv")})

    # For accumulating results across patients and subunits
    rows = []

    # Loop over every patient 
    for pt in patient_ids:

        # Construct file paths 
        file_5s = f"{project_folder}/{pt}.5S.BRD_norm_depth.tsv"
        file_45s = f"{project_folder}/{pt}.45S.BRD_norm_depth.tsv"

        # --- 5S ---
        df_5S = pd.read_csv(
            file_5s,
            sep="\t",
            header=None,
            names=["name", "pos", "cov"],
        )
        CN_5S = df_5S["cov"].mean() # Copy number estimate = average normalized depth across all positions
        rows.append({"sample_name": pt, "rDNA_gene": "5S", "CN": CN_5S}) # Add row for 5S with patient ID, gene label, and CN estimate

        # --- 45S ---
        df_45S = pd.read_csv(
            file_45s,
            sep="\t",
            header=None,
            names=["name", "pos", "cov"],
        )

        df_45S["rDNA_gene"] = annotate_45s_unit(df_45S["pos"]) # Assign subunit label per position
        df_45S = df_45S.dropna(subset=["rDNA_gene"]) # Remove positions outside the three subunit windows
        
        # Mean CN per subunit
        subunit_means = (
            df_45S.groupby("rDNA_gene", observed=True)["cov"]
            .mean() # Copy number estimate = average normalized depth across all positions
            .rename("CN")
            .reset_index()
        )

        # Iterate trough subunit means and add rows for 18S/5.8S/28S
        for _, r in subunit_means.iterrows(): 
            rows.append({"sample_name": pt, "rDNA_gene": r["rDNA_gene"], "CN": r["CN"]})

        # Total 45S = mean of the three subunit means (matches your intent)
        rows.append({"sample_name": pt, "rDNA_gene": "45S", "CN": subunit_means["CN"].mean()})

    return pd.DataFrame(rows)

def main():

    # Define command-line arguments for the script
    parser = argparse.ArgumentParser()
    parser.add_argument("--folder", required=True, help="")
    parser.add_argument("--output", required=True, help="")
    args = parser.parse_args()

    # Compute CN estimates for all patients in the folder
    CN_df = collect_rDNA_CN(args.folder)

    # Save to tsv file
    CN_df.to_csv(args.output,
                 sep="\t",
                 index=False
                 )

if __name__ == "__main__":
    main()

