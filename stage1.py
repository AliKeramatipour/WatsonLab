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


# Map flag name -> function
FLAG_FUNCTIONS = {
    '-MAINCHR': MAINCHR
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
