# Example GIFT input

This folder shows the structure and file formats the GIFT pipeline expects and
produces. **All values here are synthetic** — sample IDs, rsIDs, positions,
allele frequencies, effect sizes, and LD correlations are randomly generated.
Real ROSMAP-derived data is governed by a DUA and cannot be shared.

The example contains one region with two genes (`GENEA`, `GENEB`) and 14 SNPs
total (8 in `GENEA`, 6 in `GENEB`) across 20 fake samples. Real runs use
hundreds to thousands of SNPs per region across ~931 samples.

## Tree

```
region_example/
├── GENEA_expr.txt              # per-gene expression: sample_id <tab> value
├── GENEB_expr.txt
├── GWAS.txt                    # stacked GWAS summary stats for all SNPs in region
├── pindex.txt                  # SNP count per gene (one int per line, in gene order)
├── region_example.geno.txt     # genotype matrix: SNP info + per-sample dosages
├── snp_loc.txt                 # CHROM<tab>POS for each SNP
├── snplist.txt                 # rsIDs, one per line (matches GWAS.txt row order)
└── final_input/                # files consumed by GIFT_summary()
    ├── eQTL/
    │   ├── eQTLGENEA.txt       # eQTL summary stats for GENEA (same schema as GWAS.txt)
    │   └── eQTLGENEB.txt
    ├── GWAS.txt                # copy of region-level GWAS.txt
    ├── GWASLD.txt              # SNP x SNP LD correlation matrix from GWAS reference
    ├── LD_eQTL.txt             # SNP x SNP LD correlation matrix from eQTL reference
    ├── R_matrix.txt            # gene x gene expression correlation matrix
    ├── snplist.txt             # copy of region-level snplist.txt
    ├── pindex.txt              # copy of region-level pindex.txt
    └── output/                 # GIFT_summary() writes results here
```

## File schemas

### `GWAS.txt` and `eQTL/eQTL<GENE>.txt`

Tab-separated, with header. Columns:

| Column | Description |
|---|---|
| `chr` | Chromosome |
| `rs` | rsID |
| `ps` | Position (bp) |
| `n_mis` | Number of missing genotypes |
| `n_obs` | Sample size (GWAS: total N; eQTL: RNA-seq N) |
| `allele1` | Effect allele |
| `allelel0` | Reference allele *(note: column name uses two L's, matches GEMMA convention)* |
| `af` | Allele frequency of effect allele |
| `beta` | Effect size |
| `se` | Standard error of beta |
| `p_wald` | Wald p-value |

SNP order in `eQTL<GENE>.txt` files matches the gene's slice of `GWAS.txt`,
sized according to `pindex.txt`.

### `snplist.txt`

One rsID per line, in the same order as `GWAS.txt` rows.

### `pindex.txt`

Integer SNP count per gene, one per line, in the order genes appear in the
region. Sum of values equals the number of rows in `GWAS.txt` / `snplist.txt`.

### `snp_loc.txt`

CHROM and POS for each SNP. Tab-separated with header.

### `region_example.geno.txt`

Genotype dosage matrix. First five columns are `CHROM`, `POS`, `ID`, `REF`,
`ALT`; remaining columns are per-sample dosages (0, 1, or 2). Used to compute
LD matrices.

### `<GENE>_expr.txt`

Per-gene expression, one sample per line. Two tab-separated columns:
sample ID and expression value. No header. Sample IDs match the genotype
matrix column names.

### `R_matrix.txt`

Gene x gene expression correlation matrix. Square, symmetric, 1s on diagonal.
No header, no row names.

### `GWASLD.txt` and `LD_eQTL.txt`

SNP x SNP LD correlation matrix. Square, symmetric, 1s on diagonal. Same SNP
order as `snplist.txt`. No header, no row names. In real runs these are
~thousands x thousands; here they are 14 x 14.

## Reproducing your own example

The synthetic files here were generated with random values. To regenerate or
adapt for testing, see the project README and run the pipeline against any
test region of your choosing.
