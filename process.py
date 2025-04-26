import pandas as pd
import os
from glob import glob
from pyfaidx import Fasta
from pathlib import Path
from mygene import MyGeneInfo

# === Files ===
base_dir = '/Users/stevenmonger/Documents/Work/igor_2025/gtex/GTEx_Analysis_v10_sQTL_updated/'
fasta_path = "/Users/stevenmonger/Documents/Work/savvi/SavviDB_backend/Backend/Annotation/GRCh38.primary_assembly.genome.fa"
fasta = Fasta(fasta_path)
sqtl_files = glob(os.path.join(base_dir, "*.v10.sQTLs.signif_pairs.parquet"))
output_file = "full.tsv"

# === Functions ===
def get_acceptor_motif(chrom, strand, pos):
    if strand == '+':
        intronic = fasta[chrom][pos-18:pos+2].seq.upper()
        exonic = fasta[chrom][pos+2:pos+5].seq.upper()
    else:
        # For negative strand, reverse complement
        intronic = fasta[chrom][pos-3:pos+17].complement.reverse.seq.upper()
        exonic = fasta[chrom][pos-6:pos-3].complement.reverse.seq.upper()

    return f"{intronic}/{exonic}"

def extract_motif(row):
    try:
        return get_acceptor_motif(row['var_chr'], row['strand'], row['var_pos'])
    except KeyError:
        return "N/A"

# === Results collect ===
all_results = []
skipped_total = 0

#for file in sqtl_files:
#    tissue = os.path.basename(file).replace(".v10.sQTLs.signif_pairs.parquet", "")
#    print(f"Processing {tissue}...")
#
#    try:
#        df = pd.read_parquet(file)
#    except Exception as e:
#        print(f"Failed to load {file}: {e}")
#        continue
#
#    original_n = len(df)
#
#    # === Parse phenotype_id ===
#    split_cols = df['phenotype_id'].str.split(':', expand=True)
#    if split_cols.shape[1] != 5:
#        print(f"Skipping {tissue}, malformed phenotype_id structure.")
#        skipped_total += original_n
#        continue
#
#    df[['phen_chr', 'phen_start', 'phen_end', 'clu_with_strand', 'gene']] = split_cols
#    df['phen_start'] = df['phen_start'].astype(int)
#    df['phen_end'] = df['phen_end'].astype(int)
#
#    # === Determine strand from clu_with_strand ===
#    df['strand'] = df['clu_with_strand'].str.extract(r'clu_\d+_([-+])')
#
#    # === Parse variant_id ===
#    var_split = df['variant_id'].str.split('_', expand=True)
#    if var_split.shape[1] != 5:
#        print(f"Skipping {tissue}, malformed variant_id.")
#        skipped_total += original_n
#        continue
#
#    df[['var_chr', 'var_pos', 'ref', 'alt', 'build']] = var_split
#    df['var_pos'] = df['var_pos'].astype(int)
#
#    # === Relative position of variant to acceptor site ===
#    def compute_delta(row):
#        if row['strand'] == '+':
#            return row['phen_end'] - row['var_pos']
#        elif row['strand'] == '-':
#            return row['var_pos'] - row['phen_start']
#        else:
#            return None
#
#    df['delta'] = df.apply(compute_delta, axis=1)
#
#    # === Filter for -3 variants and C↔T transitions ===
#    df = df[df['delta'] == 3]
#    df_pos = df[df['ref'].isin(['C', 'T']) & df['alt'].isin(['C', 'T'])]
#    df_neg = df[df['ref'].isin(['A', 'G']) & df['alt'].isin(['A', 'G'])]
#    df = pd.concat([df_pos, df_neg])
#
#    # === Classify based on ref and slope direction ===
#    def classify(row):
#        slope = row['slope']
#        if pd.isna(slope):
#            return None
#        if row['ref'] == 'T' and slope > 0:
#            return 'canonical'
#        elif row['ref'] == 'T' and slope < 0:
#            return 'laggard'
#        elif row['ref'] == 'C' and slope > 0:
#            return 'laggard'
#        elif row['ref'] == 'C' and slope < 0:
#            return 'canonical'
#        if row['ref'] == 'A' and slope > 0:
#            return 'canonical'
#        elif row['ref'] == 'A' and slope < 0:
#            return 'laggard'
#        elif row['ref'] == 'G' and slope > 0:
#            return 'laggard'
#        elif row['ref'] == 'G' and slope < 0:
#            return 'canonical'
#        else:
#            return None
#
#    df['classification'] = df.apply(classify, axis=1)
#    df['tissue'] = tissue
#
#    # === Get motif ===
#    df['motif'] = df.apply(extract_motif, axis=1)
#
#    # === Save results ===
#    all_results.append(df[['variant_id', 'phenotype_id', 'ref', 'alt', 'slope', 'strand', 'classification', 'tissue', 'motif', 'af']])
#
## === Report ===
#if all_results:
#    result_df = pd.concat(all_results, ignore_index=True)
#    result_df.to_csv(output_file, sep='\t', index=False)
#    print(f"\nSaved results to {output_file}")
#
#    stats = result_df['classification'].value_counts()
#    print("\n--- Summary ---")
#    print(f"Total sQTLs at -3 with C↔T transitions: {len(result_df)}")
#    print("classification")
#    print(stats)
#else:
#    print("No variants found.")
#
#print(f"\nTotal skipped rows across all tissues: {skipped_total}")

# === Postprocessing ===
df = pd.read_csv(output_file, sep="\t")

# Create separate aggregations
phenotype_joined = df.groupby('variant_id')['phenotype_id'].apply(
    lambda x: ','.join(sorted(set(x)))
)
phenotype_count = df.groupby('variant_id')['phenotype_id'].apply(
    lambda x: len(set(x))
)
tissue_joined = df.groupby('variant_id')['tissue'].apply(
    lambda x: ','.join(sorted(set(x)))
)
tissue_count = df.groupby('variant_id')['tissue'].apply(
    lambda x: len(set(x))
)
other_columns = df.groupby('variant_id').agg({
    'strand': lambda x: ','.join(sorted(set(x))),
    'motif': lambda x: ','.join(sorted(set(x))),
    'classification': lambda x: ','.join(sorted(set(x))),
    'af': 'first'
})

def max_abs_preserving_sign(x):
    if x.empty:
        return np.nan
    abs_vals = x.abs()
    max_idx = abs_vals.idxmax()
    return x.loc[max_idx]

slope_stats = df.groupby('variant_id').agg({
    'slope': [
        max_abs_preserving_sign,
        'median'                 
    ]
})

# Flatten the multi-level column names
slope_stats.columns = ['slope_max_abs', 'slope_median']

# === Extract and validate gene_id mapping ===

df['gene_id'] = df['phenotype_id'].str.split(':').str[-1].str.split('.').str[0]

# Collapse all gene_ids per variant
gene_id_joined = df.groupby('variant_id')['gene_id'].apply(
    lambda x: ','.join(sorted(set(x)))
)

# === Query MyGene.info for unique gene IDs (after stripping version)
unique_gene_ids = sorted(set(gene_id_joined.str.split(',').sum()))
mg = MyGeneInfo()
query_results = mg.querymany(unique_gene_ids, scopes='ensembl.gene', fields='symbol', species='human')

# Build the mapping, and add manual fallback
id_to_symbol = {entry['query']: entry.get('symbol', '') for entry in query_results}
id_to_symbol['ENSG00000256427'] = 'AC010175.1'  # manually added
id_to_symbol['ENSG00000223458'] = 'LMO7DN-IT1'  # manually added
id_to_symbol['ENSG00000258414'] = 'RP11-356O9.1'  # manually added

# Map each variant's gene_id string (comma-separated) to symbols
def map_to_symbols(gene_str):
    ids = gene_str.split(',')
    symbols = [id_to_symbol.get(gid, '') for gid in ids]
    return ','.join(symbols)

gene_symbols = gene_id_joined.map(map_to_symbols)

# Combine all results
result = pd.concat([
    phenotype_joined.rename('phenotype_id'),
    phenotype_count.rename('phenotype_count'),
    tissue_joined.rename('tissue'),
    tissue_count.rename('tissue_count'),
    gene_symbols.rename('gene_symbol'),
    other_columns,
    slope_stats
], axis=1).drop(columns=['tissue','phenotype_id','strand']).sort_values(by='classification',ascending=False)

# Round float columns to 2 decimal places
float_cols = result.select_dtypes(include=['float64', 'float32']).columns
result[float_cols] = result[float_cols].round(2)

result.to_csv("aggregated.tsv", sep="\t")
