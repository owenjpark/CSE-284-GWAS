import pandas as pd
import numpy as np
from cyvcf2 import VCF
import argparse
import os

def parse_args():
    parser = argparse.ArgumentParser(
        description="Perform LD clumping on GWAS summary statistics"
    )
    parser.add_argument(
        "--clump",
        required=True,
        help="GWAS summary statistics path to use for clumping"
    )
    parser.add_argument(
        "--vcf",
        required=True,
        help="Path to input VCF file GWAS summary statistics were performed on"
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Output file path (e.g. results/clumped.tsv)"
    )
    parser.add_argument(
        "--p1",
        type=float,
        default=5e-8,
        help="P-value threshold for index SNPs (default: 5e-8)"
    )
    parser.add_argument(
        "--p2",
        type=float,
        default=0.01,
        help="P-value threshold for secondary SNPs (default: 0.01)"
    )
    parser.add_argument(
        "--kb-window",
        type=int,
        default=250,
        help="Window size in kb for clumping (default: 250)"
    )
    parser.add_argument(
        "--r2",
        type=float,
        default=0.5,
        help="LD r^2 threshold (default: 0.5)"
    )
    return parser.parse_args()

def load_gwas(input_path, p2):
    # Load results
    gwas = pd.read_csv(input_path, sep=r"\s+")
    gwas = gwas[gwas["P"] <= p2].copy()

    return gwas


def load_genotypes(vcf_path, needed_snps):
    # Load VCF and keep only needed SNPs
    vcf = VCF(vcf_path)
    
    geno_dict = {}
    seen_variants = set()
    for v in vcf:
        # Skip if not in GWAS results
        if v.ID not in needed_snps:
            continue
        
        # Check for and skip duplicates
        variant_key = (v.CHROM, v.POS, v.REF, ",".join(v.ALT))
        if variant_key in seen_variants:
            continue
        seen_variants.add(variant_key)

        # Convert to dosage
        g = v.gt_types.astype(float)
        missing_mask = (g == 2)  # UNKNOWN
        g[g == 3] = 2  # HOM_ALT -> dosage 2
        g[missing_mask] = np.nan  # UNKNOWN -> NaN

        geno_dict[v.ID] = g

    return geno_dict


def ld_r2(geno_dict, snp1, snp2):
    g1 = geno_dict[snp1]
    g2 = geno_dict[snp2]

    # Filter out NaNs and check if at least 3 overlapping samples
    mask = ~np.isnan(g1) & ~np.isnan(g2)
    if mask.sum() < 3:
        return np.nan

    # Correlation coefficient calculation
    x = g1[mask]
    y = g2[mask]

    x_mean = x.mean()
    y_mean = y.mean()

    cov = np.sum((x - x_mean) * (y - y_mean))
    var_x = np.sum((x - x_mean) ** 2)
    var_y = np.sum((y - y_mean) ** 2)

    if var_x == 0 or var_y == 0:
        return np.nan

    r = cov / np.sqrt(var_x * var_y)

    return r * r


def clump_chromosome(df_chr, geno_dict, args, lead_snps, lead_to_sp2, removed, bp_window):
    # Sort by BP for left/right scanning within chromosome
    df_chr_bp = df_chr.sort_values("BP").reset_index(drop=True)

    # Choose lead SNPs in p-value order within chromosome
    order = df_chr_bp.sort_values("P").index

    # Iterate through to clump
    for idx in order:
        row = df_chr_bp.loc[idx]
        row_snp = row["SNP"]

        # Check if already clumped
        if row_snp in removed:
            continue

        # Only SNPs passing p1 can start a clump
        if row["P"] > args.p1:
            continue

        # Add as lead SNP
        lead_snps.append(row_snp)
        lead_to_sp2[row_snp] = []
        row_bp = row["BP"]

        # Scan left side of window
        j = idx - 1
        while j >= 0:
            other = df_chr_bp.loc[j]
            other_snp = other["SNP"]

            dist_bp = row_bp - other["BP"]
            if dist_bp > bp_window:
                break

            if other_snp not in removed:
                r2 = ld_r2(geno_dict, row_snp, other_snp)
                if not np.isnan(r2) and r2 >= args.r2:
                    removed.add(other_snp)
                    lead_to_sp2[row_snp].append(other_snp)

            j -= 1
        
        # Scan right side of window
        j = idx + 1
        while j < len(df_chr_bp):
            other = df_chr_bp.loc[j]
            other_snp = other["SNP"]

            dist_bp = other["BP"] - row_bp
            if dist_bp > bp_window:
                break

            if other_snp not in removed:
                r2 = ld_r2(geno_dict, row_snp, other_snp)
                if not np.isnan(r2) and r2 >= args.r2:
                    removed.add(other_snp)
                    lead_to_sp2[row_snp].append(other_snp)

            j += 1


def clump(args):
    gwas = load_gwas(args.clump, args.p2)

    # SNPs actually needed for LD
    needed_snps = set(gwas["SNP"])

    # Create genotype dict for LD calculations
    geno_dict = load_genotypes(args.vcf, needed_snps)

    # Iterate through chromosomes to clump
    lead_snps = []
    lead_to_sp2 = {}
    removed = set()
    bp_window = args.kb_window * 1000
    for chrom, df_chr in gwas.groupby("CHR"):
        clump_chromosome(df_chr, geno_dict, args, lead_snps, lead_to_sp2, removed, bp_window)

    # Filter GWAS results to just include lead SNPs
    clumped = gwas[gwas["SNP"].isin(lead_snps)].copy()

    # Get list of other SNPs in lead SNPs group and add as column
    clumped["SP2"] = clumped["SNP"].map(lambda snp: ",".join(lead_to_sp2.get(snp, [])))

    # Save
    clumped.to_csv(args.out, sep="\t", index=False)


def main():
    args = parse_args()
    clump(args)


if __name__ == "__main__":
    main()