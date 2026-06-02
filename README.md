# BGW-TWAS-PWAS-on-AD-dementia

Bayesian Genome-wide TWAS and PWAS analyses of Alzheimer's Disease (AD)
dementia, using ROS/MAP omics data from the dorsolateral prefrontal cortex
(DLPFC) and GWAS summary statistics from Bellenguez et al. (Nature Genetics,
2022). This repository contains analysis scripts for GIFT fine-mapping,
BGW-TWAS, and BGW-PWAS.

## Background

Standard TWAS tools (PrediXcan, FUSION, TIGAR) use only *cis*-eQTL information.
Trans-eQTLs account for a substantial fraction of regulatory signal (over 30%
of genes in whole blood, and around 37% of trait-associated GWAS signals in
eQTLGen), so methods that incorporate *trans*-xQTL effects can recover risk
genes and biological pathways that cis-only methods miss.

**BGW-xWAS** (Luningham et al., AJHG, 2020) is a Bayesian variable selection
regression framework that jointly models cis- and trans- xQTL effects using
spike-and-slab priors, enabling genome-wide xWAS testing with either
individual-level or summary-level GWAS data. **GIFT** is used here for
fine-mapping of TWAS signals.

## Repository structure

```
BGW-TWAS-PWAS-on-AD-dementia/
├── README.md                  # this file
├── BGW_association/           # genome-wide TWAS/PWAS association testing (cis + trans)
│   ├── README.md
│   ├── config.example.R
│   ├── config.example.sh
│   ├── scripts/               # association_test.R, cis_association_test.R
│   ├── slurm/                 # run_association.sh, generate_geno_ss.sh
│   └── example_input/         # tiny runnable synthetic example
├── GIFT_scripts/              # GIFT input prep and helper R scripts
│   ├── config.example.R
│   ├── config.example.sh
│   ├── eQTL.R
│   ├── newGIFTGEN.R
│   ├── organization_region.R
│   └── test_input.sh
├── GIFT_pipeline/             # end-to-end GIFT fine-mapping pipeline
│   ├── README.md
│   ├── config.example.R
│   ├── runall.sh
│   └── xWASGIFT_pipeline.R
└── GIFT_example_input/        # example region inputs for the GIFT pipeline
    ├── README.md
    └── regionexample/
        ├── final_input/
        │   ├── eQTL/           # eQTLGENEA.txt, eQTLGENEB.txt
        │   ├── GWAS.txt, GWASLD.txt, LD_eQTL.txt, R.txt
        │   └── snplist.txt, snploc.txt
        ├── GENEA_expr.txt, GENEB_expr.txt
        ├── GWAS.txt, pindex.txt, region_example.geno.txt
        └── snplist.txt, snploc.txt
```

Each component holds its own `config.example.*` templates. Copy them to
`config.R` / `config.sh` and edit the paths for your environment; the copies
are gitignored so no absolute, DUA-governed paths are committed. See the README
inside each component for setup and run instructions (start with
`BGW_association/README.md` for the association testing step).

## Data

- **Transcriptomics:** ROS/MAP DLPFC bulk RNA-seq (n = 931)
- **Proteomics:** ROS/MAP DLPFC TMT mass spectrometry (n = 716)
- **Genotypes:** ROS/MAP whole genome sequencing
- **GWAS:** Bellenguez et al. 2022 (111,326 AD cases / 677,663 controls)

Data are accessed under the relevant ROS/MAP and consortia data use agreements
and are **not** included in this repository. Example inputs provided here are
small synthetic files for format illustration and testing only.

## Status

Active development. Results, figures, and final documentation will be added in
later commits.

## Citation

Luningham JM, Chen J, Tang S, De Jager PL, Bennett DA, Buchman AS, Yang J.
*Bayesian Genome-wide TWAS Method to Leverage both cis- and trans-eQTL
Information through Summary Statistics.* Am J Hum Genet. 2020;107(4):714-726.
[doi:10.1016/j.ajhg.2020.08.022](https://doi.org/10.1016/j.ajhg.2020.08.022)

## Contact

Arpit Ramani, Yang Lab, Emory University School of Medicine.
