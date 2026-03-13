from cyvcf2 import VCF
import numpy as np
import pandas as pd
import argparse
import os


def parse_args():
    parser = argparse.ArgumentParser(description="Compute PCA from a VCF using the GRM.")
    parser.add_argument(
        "--num_pcs",
        type=int,
        default=None,
        help="Number of principal components to compute (default: all)"
    )
    parser.add_argument(
        "--vcf",
        required=True,
        help="Path to input VCF file"
    )
    parser.add_argument(
        "--outdir",
        default=".",
        help="Output directory (default: current directory)"
    )
    return parser.parse_args()


def load_genotype_matrix(vcf_path):
    # Read VCF
    vcf = VCF(vcf_path)
    samples = vcf.samples

    # Collect SNPs into list
    genotype_list = []
    for variant in vcf:
        g = np.array([gt[0] + gt[1] if gt[0] >= 0 and gt[1] >= 0 else np.nan for gt in variant.genotypes])
        genotype_list.append(g)

    # Shape: (n_samples, n_snps)
    X = np.vstack(genotype_list).T
    return samples, X


def compute_pca_from_genotypes(X, num_pcs=None):
    n_samples = X.shape[0]

    # If num_pcs not provided, return all PCs
    k = num_pcs if num_pcs is not None else n_samples

    # Check if not enough PCs
    if k < 1 or k > n_samples:
        raise ValueError(f"num_pcs must be between 1 and {n_samples}, got {k}")

    # Remove SNPs with no data
    valid = ~np.isnan(X).all(axis=0)
    X = X[:, valid]

    # Estimate allele frequencies
    p = np.nanmean(X, axis=0) / 2

    # Expected variance under Hardy-Weinberg
    var = 2 * p * (1 - p)

    # Remove monomorphic SNPs
    valid = (~np.isnan(var)) & (var > 0)
    X = X[:, valid]
    p = p[valid]
    var = var[valid]

    # Standardize genotypes
    Z = (X - 2 * p) / np.sqrt(var)

    # Impute missing values with mean (0 after standardization)
    Z = np.where(np.isnan(Z), 0, Z)

    # Compute GRM
    GRM = (Z @ Z.T) / Z.shape[1]

    # PCA via eigendecomposition
    eigenvals_all, eigenvecs_all = np.linalg.eigh(GRM)

    # Take top k PCs (descending)
    eigenvals = eigenvals_all[-k:][::-1]
    PCs = eigenvecs_all[:, -k:][:, ::-1]

    return eigenvals, PCs


def save_results(samples, pcs, eigenvals, outdir):
    k = pcs.shape[1]

    # Ensure output directory exists
    os.makedirs(outdir, exist_ok=True)

    # Save eigenvectors
    pcs_df = pd.DataFrame(pcs, columns=[f"PC{i+1}" for i in range(k)])
    pcs_df.insert(0, "IID", samples)
    pcs_df.insert(0, "FID", samples)
    pcs_df.to_csv(os.path.join(outdir, "eigenvec.txt"), sep=" ", header=False, index=False)

    # Save eigenvalues
    eval_df = pd.DataFrame({"PC": [f"PC{i+1}" for i in range(k)],"eigenvalue": eigenvals})
    eval_df.to_csv(os.path.join(outdir, "eigenval.txt"), sep="\t", index=False)


def main():
    args = parse_args()
    samples, X = load_genotype_matrix(args.vcf)
    eigenvals, pcs = compute_pca_from_genotypes(X, args.num_pcs)
    save_results(samples, pcs, eigenvals, args.outdir)


if __name__ == "__main__":
    main()