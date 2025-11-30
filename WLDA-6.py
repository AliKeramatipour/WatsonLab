import sys
import os
import pandas as pd

# Keep only these chromosomes
MAIN_CHROMS = {
    'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
    'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
    'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrMT'
}

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 WLDA-6.py <input_file>")
        sys.exit(1)

    input_file = sys.argv[1]

    if not os.path.exists(input_file):
        print(f"Error: input file '{input_file}' not found")
        sys.exit(1)

    # Prepare output paths
    input_dir = os.path.dirname(os.path.abspath(input_file))
    output_dir = os.path.join(input_dir, "WLDA-6")
    os.makedirs(output_dir, exist_ok=True)

    output_file = os.path.join(output_dir, "output.txt")

    # Load file
    df = pd.read_csv(input_file, sep="\t", dtype=str)

    if "Chr" not in df.columns:
        print("Error: 'Chr' column not found in the file.")
        sys.exit(1)

    # Keep only rows where first column (Chr) is a main chromosome
    filtered = df[df["Chr"].isin(MAIN_CHROMS)]

    # Write output
    filtered.to_csv(output_file, sep="\t", index=False)

    removed = len(df) - len(filtered)

    print("✓ Finished")
    print(f"Input: {input_file}")
    print(f"Output folder: {output_dir}")
    print(f"Output file: {output_file}")
    print(f"Rows removed: {removed}")

if __name__ == "__main__":
    main()