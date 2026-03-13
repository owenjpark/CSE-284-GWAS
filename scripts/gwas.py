from cyvcf2 import VCF
import pandas as pd
import numpy as np
from scipy import stats
import argparse
import os


def parse_args():
    parser = argparse.ArgumentParser(description="Run GWAS with optional covariates.")
    parser.add_argument(
        "--covar",
        action="store_true",
        help="Include covariates (PCs) in the regression."
    )
    parser.add_argument(
        "--vcf",
        required=True,
        help="Path to input VCF file"
    )
    parser.add_argument(
        "--phen",
        required=True,
        help="Path to phenotype file."
    )
    parser.add_argument(
        "--pcs",
        help="Path to eigenvec/PC file."
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Output file path (e.g. results/gwas.tsv)"
    )
    parser.add_argument(
        "--maf",
        type=float,
        default=0.01,
        help="Minor allele frequency threshold (default: 0.01)"
    )
    return parser.parse_args()


def load_phenotypes(vcf, phen_path):
    # Reorder phenotypes samples to match order of VCF samples
    phen = pd.read_csv(phen_path, sep="\t", header=None)
    phen.columns = ["FID", "IID", "PHEN"]
    phen = phen.set_index("IID")
    samples = vcf.samples
    y = phen.loc[samples]["PHEN"].values.astype(float)

    return samples, y


def load_covariates(samples, pcs_path):
    # Reorder PCs to match order of VCF samples
    pcs = pd.read_csv(pcs_path, sep=r"\s+", header=None)

    pcs.columns = ["FID", "IID"] + [f"PC{i}" for i in range(1, pcs.shape[1] - 1)]
    pcs = pcs.set_index("IID")
    pc_matrix = pcs.loc[samples].drop(columns=["FID"]).values.astype(float)

    # Add intercept/bias term (columns of 1s added to left side)
    intercept = np.ones(len(pc_matrix))

    C = np.column_stack([intercept, pc_matrix])

    # Creates a boolean mask indicating which rows have no missing values in C precomputed for later use: (n_samples,)
    covar_complete = ~np.any(np.isnan(C), axis=1)

    return C, covar_complete


def genotype_to_dosage(v):
    # Convert genotype to dosage-like numeric
    g = v.gt_types.astype(float)
    missing_mask = (g == 2) # UNKNOWN
    g[g == 3] = 2 # HOM_ALT -> dosage 2
    g[missing_mask] = np.nan # UNKNOWN -> NaN

    return g


def run_gwas(args):
    # Read VCF and phenotype data
    vcf = VCF(args.vcf)
    samples, y = load_phenotypes(vcf, args.phen)

    # Optional covariate setup
    C = None
    covar_complete = None
    y_resid = None
    if args.covar:
        if args.pcs is None:
            raise ValueError("--pcs is required when using --covar")
        C, covar_complete = load_covariates(samples, args.pcs)
        covar_coefs_y = np.linalg.solve(C.T @ C, C.T @ y)
        y_resid = y - C @ covar_coefs_y

    # Set output filename from prefix
    prefix = args.out
    out_path = f"{prefix}.assoc.linear"

    # Ensure output directory exists
    out_dir = os.path.dirname(out_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    # Set up output file
    with open(out_path, "w") as out:
        out.write("CHR\tBP\tSNP\tBETA\tP\tN\n")

        # Iterate through each SNP
        for v in vcf:
            # Convert genotype to dosage-like numeric
            g = genotype_to_dosage(v)

            # Masking behavior based on if using covar
            if args.covar:
                mask = ~np.isnan(g) & ~np.isnan(y_resid) & covar_complete
                g2 = g[mask]
                y2 = y_resid[mask]
                C2 = C[mask]
                if len(g2) < C2.shape[1] + 2:
                    continue
            else:
                mask = ~np.isnan(g) & ~np.isnan(y)
                g2 = g[mask]
                y2 = y[mask]
                if len(g2) < 3:
                    continue

            # Calculate MAF aka p
            p = g2.mean() / 2

            # Flip allele to match PLINK using minor allele, not ALT allele as effect allele
            if p > 0.5:
                g2 = 2 - g2
                p = 1 - p

            # Don't include SNP if MAF less than specified value
            if p < args.maf:
                continue
            
            # Regression behavior based on if using covar
            if args.covar:
                # Compute g-residuals (removes the effect of the PCs and intercept from the genotypes)
                covar_coefs_g = np.linalg.solve(C2.T @ C2, C2.T @ g2)
                g_used = g2 - C2 @ covar_coefs_g
                y_used = y2

                # Regression on filtered samples
                denom = g_used @ g_used
                if denom == 0: # Skip monomorphic SNPs (only 1 allele present in all samples)
                    continue

                beta = (g_used @ y_used) / denom
                dof = len(y_used) - C2.shape[1] - 1

                resid = y_used - beta * g_used
                se = np.sqrt((resid @ resid) / dof / denom)
            else:
                y0 = y2 - y2.mean()
                g0 = g2 - g2.mean()

                denom = g0 @ g0
                if denom == 0: # Skip monomorphic SNPs (only 1 allele present in all samples)
                    continue

                beta = (g0 @ y0) / denom
                
                dof = len(y2) - 2
                resid = y0 - beta * g0
                se = np.sqrt((resid @ resid) / dof / denom)

            t = beta / se
            pval = 2 * stats.t.sf(abs(t), dof)

            # Write to file
            out.write(f"{v.CHROM}\t{v.POS}\t{v.ID}\t{beta}\t{pval}\t{len(g2)}\n")


def main():
    args = parse_args()
    run_gwas(args)


if __name__ == "__main__":
    main()