# BGW association test

Genome-wide TWAS / PWAS association testing using BGW-xWAS weights (cis and
trans xQTL effects) against GWAS summary statistics. Given per-gene BGW-xWAS
weight files, genotype dosages, and GWAS summary stats, the scripts compute a
burden-style Z statistic and a TWAS/PWAS p-value for each gene, run per
chromosome as a SLURM array job.

This module is self-contained and path-portable. All site-specific paths live
in `config.R` / `config.sh`, which are gitignored, so the scripts themselves
contain no absolute paths.

## Layout

```
BGW_association/
├── config.example.R       # R-side data paths (copy to config.R)
├── config.example.sh      # shell/SLURM paths and compute settings (copy to config.sh)
├── scripts/
│   ├── association_test.R       # full test (cis + trans weights)
│   └── cis_association_test.R   # cis-only test (Trans == 0), writes zero-cis gene list
├── slurm/
│   ├── run_association.sh        # array job: run the association test per chromosome
│   └── generate_geno_ss.sh       # array job: build .geno.txt / .SS.txt via Estep_mcmc
└── example_input/         # tiny synthetic example (one gene, runnable end to end)
    ├── Gene_lists/Gene_list_1.txt
    ├── weights/EXAMPLEGENE/EXAMPLEGENE_BGW_xQTL_weights.txt
    ├── Test_Geno/EXAMPLEGENE.geno.txt
    ├── Test_GWAS_SS/EXAMPLEGENE.SS.txt
    └── Results/
```

## Setup

1. Copy and edit the two config templates:

   ```
   cp config.example.R  config.R
   cp config.example.sh config.sh
   ```

   Set the data directories in `config.R` (`WEIGHT_DIR`, `GENO_DIR`, `SS_DIR`,
   `RESULT_DIR`) and the matching paths plus compute settings in `config.sh`.
   Keep `GENO_DIR`/`SS_DIR` in `config.R` consistent with `GENO_OUT_DIR`/`SS_OUT_DIR`
   in `config.sh`.

2. Set `PROJECT_DIR` in `config.sh` to the absolute path of this module so the
   R scripts can locate `config.R` via the exported `CONFIG_R` variable.

## Running

Per-chromosome array jobs (defaults to chromosomes 1 to 22). Submit from the
module root so the `logs/` paths resolve, and edit the `#SBATCH --partition`
line for your cluster:

```
mkdir -p logs

# Optional: build genotype dosage and summary-stat files first
sbatch slurm/generate_geno_ss.sh

# Run the association test
sbatch slurm/run_association.sh
```

To run the cis-only test, edit `slurm/run_association.sh` and switch the
`R_SCRIPT` line to `cis_association_test.R`.

You can also call a script directly (outside SLURM):

```
CONFIG_R=$PWD/config.R Rscript scripts/association_test.R <gene_list_file> <chr> <num_cores>
```

## Trying the example

Point a throwaway `config.R` at `example_input/` and run chromosome 1:

```
cat > /tmp/config.R <<'EOF'
WEIGHT_DIR <- "example_input/weights"
GENO_DIR   <- "example_input/Test_Geno"
SS_DIR     <- "example_input/Test_GWAS_SS"
RESULT_DIR <- "example_input/Results"
EOF

CONFIG_R=/tmp/config.R Rscript scripts/association_test.R \
  example_input/Gene_lists/Gene_list_1.txt 1 1
```

This writes `example_input/Results/CHR_1_assoc.tsv` with a single gene row
(EXAMPLEGENE), producing a finite ZStat and p-value.

## Input formats

Weight file (`<gene>_BGW_xQTL_weights.txt`), 13 columns, header optional:
`CHROM POS ID REF ALT MAF Trans CPP Beta mBeta ChisqTest Pval_svt Rank`
(`Trans`: 0 = cis, 1 = trans).

Genotype dosage (`<gene>.geno.txt`): `CHROM POS ID REF ALT` followed by one
column per sample (dosage values).

GWAS summary stats (`<gene>.SS.txt`):
`CHROM POS rsID REF ALT Beta_gwas Beta_sd_gwas PVAL_gwas`.

SNPs are matched across files by a `CHROM:POS:REF:ALT` ID (a `chr` prefix on
chromosome names is stripped automatically), with allele-flip handling on the
GWAS side.

## Output

`CHR_<chr>_assoc.tsv` (full) or `CHR_<chr>_cis_assoc.tsv` (cis-only), with
columns: `GeneName, ZStat, PVAL_twas, N_with_weight, Ncis_with_weight,
N_tested, Ncis_tested, sum_cpp_cis, sum_cpp_trans`. The cis-only script also
writes `CHR_<chr>_zero_cis_genes.txt` listing genes whose weights had no cis
SNPs (these are dropped from the results table).

## Requirements

R with `data.table`, `tidyverse`, `foreach`, and `doParallel`. The genotype
generation step additionally needs the BGW-xWAS `Estep_mcmc` binary.
