# GIFT pipeline (xWAS)

A generic 5-step pipeline that prepares GIFT inputs from BGW xWAS association
results and submits per-region GIFT fine-mapping jobs. Works for any xWAS
analysis (TWAS, PWAS, etc.) — only the config values change.

## Order of operations

| Step | Purpose | Inputs | Outputs |
|------|---------|--------|---------|
| 1 | Merge per-chromosome association files, apply genomic control adjustment, select significant genes | `ASSOC_DIR/<prefix><chr><suffix>` | `merged_gc_adjusted.tsv`, `merged_gc_adjusted_significant.tsv` |
| 2 | Match significant genes to BGW xQTL weights and per-gene summary stats | `merged_gc_adjusted_significant.tsv`, `WEIGHTS_BASE/`, `SS_DIRS` | `matched_weights/<gene>_matched.txt`, `missing_genes_in_matched.txt` |
| 3 | Combine all matched per-gene files into one table | `matched_weights/*_matched.txt` | `ALL_GENES_matched_combined_significant.txt` |
| 4 | Define LD-merged genomic regions (genes within 1Mb collapse into one region) | `ALL_GENES_matched_combined_significant.txt`, `EXPR_FILE` | `INPUT_BASE/region_from_matched.txt` |
| 5 | Build per-region GIFT inputs (GWAS, eQTL, LD, R matrix) and write SLURM submission scripts | region file, weights, GEMMA, geno files | `INPUT_BASE/region<N>/`, `INPUT_BASE/GIFT_running/` |

## Running the pipeline

### Run everything in one job

```bash
sbatch run_all.sh
```

This runs steps 1–5 sequentially. Step 5 generates per-region SLURM scripts but
does not submit them.

### Run a single step or a range

```bash
bash run_all.sh 1 1     # just step 1
bash run_all.sh 3 5     # steps 3 through 5
```

### Run a step directly

```bash
Rscript xwas_gift_pipeline.R 2
```

Bypasses the SLURM wrapper. Useful for debugging.

## Submitting the per-region GIFT jobs

Step 5 generates a submission script at `$INPUT_BASE/GIFT_running/submit_all_regions.sh`.
After the pipeline completes:

```bash
bash $INPUT_BASE/GIFT_running/submit_all_regions.sh
```

## Configuration

All paths and settings are read from `config.R` (R-side) and `config.sh`
(shell-side) in the repo root. See `config.example.R` for the full list of
variables and what each one means.

Key settings:

- `CIS_TRANS_FILTER` — `"cis"`, `"trans"`, or `""` (no filter)
- `P_CUTOFF` — significance threshold for selecting genes
- `DO_PRUNE` / `PRUNE_REGION_IDS` — optional LD pruning
- `N_SAMPLES` — sample size of the xQTL reference (e.g. RNA-seq or proteomics N)

## Notes

- Step 1 applies a genomic control correction (λ = median χ² / qchisq(0.5, 1))
  before significance filtering.
- Step 2 detects column names in summary stat files via regex, so it tolerates
  variations like `CHR`/`chr`/`CHROM`, `POS`/`ps`/`BP`, etc.
- Step 5 handles allele flipping between weights and reference GWAS files,
  removes monomorphic SNPs, and prunes SNPs missing from the genotype matrix.
- `GWASLD.txt` and `LD_eQTL.txt` are computed from the same genotype matrix in
  this pipeline. GIFT treats them as conceptually distinct LD references; here
  the same source is used for both.
