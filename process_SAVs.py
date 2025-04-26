import pandas as pd
import gzip

variants = pd.read_csv("savvi_variants_feb_2025_pos_only.csv")  
variants = variants.rename(columns={"start": "pos"})

# === Parse GENCODE GTF to find 3' splice sites ===
def parse_gtf_acceptors(gtf_file):
    acceptors = []
    with gzip.open(gtf_file, "rt") if gtf_file.endswith(".gz") else open(gtf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split('\t')
            if fields[2] != 'exon':
                continue
            chrom, start, end, strand = fields[0], int(fields[3]), int(fields[4]), fields[6]
            attrs = {k.strip(): v.strip('"') for part in fields[8].split(';') if part for k, v in [part.strip().split(' ', 1)]}
            gene = attrs.get("gene_name") or attrs.get("gene_id")
            acceptor_pos = start if strand == '+' else end
            acceptors.append({"chrom": chrom, "acceptor_pos": acceptor_pos, "strand": strand, "gene": gene})
    return pd.DataFrame(acceptors)

# === Classify variants as laggard/canonical ===
def classify_variants(variants, acceptors):
    merged = variants.merge(acceptors, on="gene", how="inner")

    def is_minus3(row):
        return (row['pos'] == row['acceptor_pos'] - 3) if row['strand'] == '+' else (row['pos'] == row['acceptor_pos'] + 3)

    minus3 = merged[merged.apply(is_minus3, axis=1)].copy()

    def classify(row):
        ref, alt = row['ref'].upper(), row['alt'].upper()
        if row['strand'] == '+':
            if ref == 'T' and alt == 'C':
                return 'laggard'
            elif ref == 'C' and alt == 'T':
                return 'canonical'
        else:
            # Convert to reverse complement to interpret as if +
            rc = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
            ref_rc = rc.get(ref, ref)
            alt_rc = rc.get(alt, alt)
            if ref_rc == 'T' and alt_rc == 'C':
                return 'laggard'
            elif ref_rc == 'C' and alt_rc == 'T':
                return 'canonical'
        return None  # Not C↔T transition

    minus3['classification'] = minus3.apply(classify, axis=1)
    minus3 = minus3[minus3['classification'].notnull()]  # Remove non-transitions

    return minus3[['chrom', 'pos', 'ref', 'alt', 'gene', 'strand', 'acceptor_pos', 'classification']].drop_duplicates(keep="first")

# === Run ===
gtf_path = "gencode.v47.primary_assembly.basic.annotation.gtf.gz"  # Adjust as needed
acceptors_df = parse_gtf_acceptors(gtf_path)
classified_variants = classify_variants(variants, acceptors_df)
classified_variants.to_csv("classified_splice_variants.csv", index=False)

# === Print stats ===
print("\n✅ Classification Summary:")
print(f"Total classified variants at -3 position: {len(classified_variants)}")
for label in ['laggard', 'canonical']:
    count = (classified_variants['classification'] == label).sum()
    unique = classified_variants[classified_variants['classification'] == label][['chrom', 'pos']].drop_duplicates().shape[0]
    print(f"  {label.capitalize()}: {count} total, {unique} unique positions")

