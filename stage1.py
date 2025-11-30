#!/usr/bin/env python3
import sys
import os
import pandas as pd

def MAINCHR(df): # WLDA-6 Keep only main chromosomes
    MAIN_CHROMS = {
        'chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',
        'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',
        'chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrMT'
    }
    if 'Chr' not in df.columns:
        print("Error: 'Chr' column not found for WLDA6.")
        sys.exit(1)
    return df[df['Chr'].isin(MAIN_CHROMS)].reset_index(drop=True)

def SPLITGENE(df): # WLDA-7
    df["Gene.refGene"] = df["Gene.refGene"].str.split(",")
    df = df.explode("Gene.refGene")
    return df.reset_index(drop=True)

def GNOMAD0(df): # WLDA-9
    df["gnomad41_genome_AF"] = df["gnomad41_genome_AF"].replace(".", "0")
    df["gnomad41_exome_AF"]  = df["gnomad41_exome_AF"].replace(".", "0")
    return df

def ACMG(df):# WLDA-10
    # Find where InterVar_automated is
    idx = df.columns.get_loc("InterVar_automated")

    # Collect the 28 ACMG evidence column names (the next 28 columns)
    acmg_cols = df.columns[idx+1 : idx+29]

    # Create empty ACMG column initially
    df.insert(idx + 1, "ACMG", "")

    # Build ACMG string row-by-row
    def summarize(row):
        names = []
        for col in acmg_cols:
            val = row[col]
            if val == "1":
                names.append(col)
        return ";".join(names) if names else ""

    df["ACMG"] = df.apply(summarize, axis=1)

    # Remove the 28 columns afterwards
    df = df.drop(columns=acmg_cols)

    return df

def CHR_POS_REF_ALT(df): # WLDA-11
    cols_to_drop = ["Chr", "Start", "End", "Ref", "Alt"]
    df = df.drop(columns=cols_to_drop, errors="ignore")

    needed = ["Otherinfo4", "Otherinfo5", "Otherinfo7", "Otherinfo8"]

    for col in needed:
        if col not in df.columns:
            print(f"Error: required column '{col}' not found.")
            sys.exit(1)

    new_front = df[needed].copy()

    rename_map = {
        "Otherinfo4": "Chr",
        "Otherinfo5": "position",
        "Otherinfo7": "Reference",
        "Otherinfo8": "Alternate"
    }
    new_front = new_front.rename(columns=rename_map)
    df = df.drop(columns=needed)
    df = pd.concat([new_front, df], axis=1)
    
    return df

def ZYGO(df): # WLDA-12
    def parse(x):
        parts = x.split(':')
        ref, alt = map(int, parts[1].split(','))
        cov = int(parts[2])
        vaf = alt / cov if cov > 0 else 0

        if vaf < 0.25:
            zy = "FP/HET"
        elif vaf < 0.75:
            zy = "HET"
        elif vaf < 0.85:
            zy = "HET/HOM"
        else:
            zy = "HOM"

        return pd.Series([zy, cov, ref, alt, vaf])

    df[["Zygosity", "Coverage", "RefReads", "AltReads", "VAF"]] = (
        df["Otherinfo13"].apply(parse)
    )

    # Insert after Alternate
    alt_idx = df.columns.get_loc("Alternate") + 1
    cols_new = ["Zygosity", "Coverage", "RefReads", "AltReads", "VAF"]
    df = df[df.columns[:alt_idx].tolist() + cols_new + df.columns[alt_idx:].tolist()]

    # Drop unwanted
    drop_cols = [
        "Otherinfo1","Otherinfo2","Otherinfo3","Otherinfo9",
        "Otherinfo10","Otherinfo11","Otherinfo12","Otherinfo13"
    ]
    return df.drop(columns=[c for c in drop_cols if c in df])

def HGVSC_P(df): # WLDA-13
    def merge_values(row):
        ref = row["AAChange.refGene"]
        ens = row["AAChange.ensGene"]

        if ref == "." and ens == ".":
            return "."
        if ref == ".":
            return ens
        if ens == ".":
            return ref
        return f"{ref} |{ens}"

    df["HgvsC&HgvsP"] = df.apply(merge_values, axis=1)

    # Remove old columns and reorder so the new column replaces refGene
    col_pos = df.columns.get_loc("AAChange.refGene")
    df.drop(columns=["AAChange.refGene", "AAChange.ensGene"], inplace=True)

    cols = list(df.columns)
    # Insert new col at the old position
    cols.insert(col_pos, cols.pop(cols.index("HgvsC&HgvsP")))
    df = df[cols]

    return df

def TRANSCRIPT(df): # WLDA-14
    def merge_values(row):
        ref = row["GeneDetail.refGene"]
        ens = row["GeneDetail.ensGene"]

        if ref == "." and ens == ".":
            return "."
        if ref == ".":
            return ens
        if ens == ".":
            return ref
        return f"{ref} |{ens}"

    # Create new merged column
    df["Transcripts"] = df.apply(merge_values, axis=1)

    # Insert new column at old refGene position
    pos = df.columns.get_loc("GeneDetail.refGene")

    df.drop(columns=["GeneDetail.refGene", "GeneDetail.ensGene"], inplace=True)

    cols = list(df.columns)
    cols.insert(pos, cols.pop(cols.index("Transcripts")))
    df = df[cols]

    return df

def DELCOL(df): # WLDA-15
    cols_to_delete = [
        "ONCDN", "ONCDISDB", "ONCREVSTAT", "ONC",
        "SCIDN", "SCIDISDB", "SCIREVSTAT", "SCI",

        "gnomad41_exome_fafmax_faf95_max",
        "gnomad41_exome_fafmax_faf99_max",
        "gnomad41_exome_AF_afr",
        "gnomad41_exome_AF_amr",
        "gnomad41_exome_AF_asj",
        "gnomad41_exome_AF_eas",
        "gnomad41_exome_AF_fin",
        "gnomad41_exome_AF_mid",
        "gnomad41_exome_AF_nfe",
        "gnomad41_exome_AF_remaining",
        "gnomad41_exome_AF_sas",

        "REGENERON_ALL_AF",
        "REGENERON_ALL_AC",
        "REGENERON_ALL_AN",

        "MCAP",
        "REVEL",
        "CLNALLELEID",

        "gnomad41_genome_fafmax_faf95_max",
        "gnomad41_genome_fafmax_faf99_max",
        "gnomad41_genome_AF_afr",
        "gnomad41_genome_AF_ami",
        "gnomad41_genome_AF_amr",
        "gnomad41_genome_AF_asj",
        "gnomad41_genome_AF_eas",
        "gnomad41_genome_AF_fin",
        "gnomad41_genome_AF_mid",
        "gnomad41_genome_AF_nfe",
        "gnomad41_genome_AF_remaining",
        "gnomad41_genome_AF_sas",

        "Func.ensGene",
        "Gene.ensGene",
        "ExonicFunc.ensGene",
        "CLNREVSTAT"
    ]

    existing = [c for c in cols_to_delete if c in df.columns]
    df = df.drop(columns=existing)

    return df

def RENAME(df): # WLDA-16
    rename_map = {
        "Gene.refGene": "Gene Name",
        "Func.refGene": "Genomic_Context",
        "ExonicFunc.refGene": "Consequence",
    }

    # Only rename columns that actually exist
    valid_map = {old: new for old, new in rename_map.items() if old in df.columns}

    df = df.rename(columns=valid_map)
    return df

def REMOVECHR(df): # WLDA-17
    if "Chr" not in df.columns:
        return df   # Nothing to do

    def strip_chr(val):
        if isinstance(val, str) and val.startswith("chr"):
            return val[3:]  # remove 'chr'
        return val

    df["Chr"] = df["Chr"].apply(strip_chr)

    # Normalise special cases
    df["Chr"] = df["Chr"].replace({"MT": "M", "m": "M", "Mt": "M"})

    return df

def REORDER(df): # WLDA-18
    final_order = [
        "Chr",
        "Position",
        "Reference",
        "Alternate",
        "Zygosity",
        "Coverage",
        "AltReads",
        "RefReads",
        "Vaf",
        "Gene Name",
        "Genomic_Context",
        "Consequence",
        "HgvsC&HgvsP",
        "Transcripts",
        "InterVar_automated",
        "ACMG",
        "CLNSIG",
        "CLNDN",
        "CLNDISDB",
        "gnomad41_genome_AF",
        "gnomad41_genome_AF_raw",
        "gnomad41_genome_AF_XX",
        "gnomad41_genome_AF_XY",
        "gnomad41_genome_AF_grpmax",
        "gnomad41_genome_faf95",
        "gnomad41_genome_faf99",
        "gnomad41_exome_AF",
        "gnomad41_exome_AF_raw",
        "gnomad41_exome_AF_XX",
        "gnomad41_exome_AF_XY",
        "gnomad41_exome_AF_grpmax",
        "gnomad41_exome_faf95",
        "gnomad41_exome_faf99",
        "SIFT_score",
        "SIFT_pred",
        "Polyphen2_HDIV_score",
        "Polyphen2_HDIV_pred",
        "Polyphen2_HVAR_score",
        "Polyphen2_HVAR_pred",
        "LRT_score",
        "LRT_pred",
        "MutationTaster_score",
        "MutationTaster_pred",
        "MutationAssessor_score",
        "MutationAssessor_pred",
        "FATHMM_score",
        "FATHMM_pred",
        "RadialSVM_score",
        "RadialSVM_pred",
        "LR_score",
        "LR_pred",
        "VEST3_score",
        "CADD_raw",
        "CADD_phred",
        "GERP++_RS",
        "phyloP46way_placental",
        "phyloP100way_vertebrate",
        "SiPhy_29way_logOdds",
    ]

    # Check for missing columns
    missing = [col for col in final_order if col not in df.columns]

    # Check for extra columns
    extra = [col for col in df.columns if col not in final_order]

    if extra:
        raise ValueError(
            "\n❌ ERROR: Unexpected columns exist in the file.\n"
            f"❗Extra columns found:\n{extra}\n"
            "The workflow has been stopped."
        )

    # Optional: warn if expected columns are missing (doesn't stop)
    if missing:
        print("\n WARNING: Some expected columns are missing:")
        print(missing)

    # Reorder only the columns that exist (missing ones are ignored)
    cols_to_use = [col for col in final_order if col in df.columns]

    return df[cols_to_use]

# Map flag name -> function
FLAG_FUNCTIONS = {
    '-MAINCHR': MAINCHR,
    '-SPLITGENE': SPLITGENE,
    '-GNOMAD0' : GNOMAD0,
    '-ACMG': ACMG,
    '-CHR-POS-REF-ALT': CHR_POS_REF_ALT,
    '-ZYGO': ZYGO,
    '-HGVSC_P': HGVSC_P,
    '-TRANSCRIPT': TRANSCRIPT,
    '-DELCOL' : DELCOL,
    '-RENAME' : RENAME,
    '-REMOVECHR': REMOVECHR,
    '-REORDER': REORDER
}

def main():
    if len(sys.argv) < 3:
        print("Usage: python3 stage1.py <input_file> <flags...>")
        sys.exit(1)

    input_file = sys.argv[1]
    flags = sys.argv[2:]

    if not os.path.exists(input_file):
        print(f"Error: input file '{input_file}' not found")
        sys.exit(1)
    else:
        print(f"Processing: {input_file}")

    for flag in flags:
        if flag not in FLAG_FUNCTIONS:
            print(f"Warning: unknown flag '{flag}', skipping")
            sys.exit(1)

    # Load once
    df = pd.read_csv(input_file, sep='\t', dtype=str)

    # Apply flags in order
    for flag in flags:
        df = FLAG_FUNCTIONS[flag](df)
        print(f"Applied: {flag}")

    # Prepare output folder
    input_dir = os.path.dirname(os.path.abspath(input_file))
    output_dir = os.path.join(input_dir, 'stage1')
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, 'output.txt')

    # Save final output
    df.to_csv(output_file, sep='\t', index=False)
    print("✓ Finished")
    print(f"Applied flags: {flags}")

if __name__ == "__main__":
    main()
