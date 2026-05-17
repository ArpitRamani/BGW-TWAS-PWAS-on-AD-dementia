if (file.exists("config.R")) {
  source("config.R")
} else if (file.exists("../config.R")) {
  source("../config.R")
} else {
  stop("config.R not found. Copy config.example.R to config.R and edit paths.")
}

options(scipen = 999)

expr    <- read.delim(EXPR_FILE)
regions <- read.delim(REGION_FILE)

for (i in 1:length(unique(regions[, "region"]))) {
  region <- regions[which(regions[, "region"] == i), ]
  region <- region[order(region$GeneName), ]
  for (j in 1:nrow(region)) {
    gene_expr <- expr[which(expr[, "GeneName"] == region[j, "GeneName"]), ]
    gene_expr <- t(gene_expr[, 6:ncol(gene_expr)])
    gene_expr <- cbind(rownames(gene_expr), gene_expr)
    gene_expr[, 1] <- gsub("\\.", "-", gene_expr[, 1])
    gene_expr <- na.omit(gene_expr)
    write.table(
      gene_expr,
      file.path(INPUT_BASE, paste0("region", i), paste0(region[j, 1], "_expr.txt")),
      sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE
    )
  }
}

options(scipen = 999)
eQTL <- c()
regions <- read.delim(REGION_FILE)

for (i in 1:length(unique(regions[, "region"]))) {
  region <- regions[which(regions[, "region"] == i), ]
  region <- region[order(region$GeneName), ]
  for (j in 1:nrow(region)) {
    region_dir <- file.path(INPUT_BASE, paste0("region", i))
    gene_dir   <- file.path(region_dir, region[j, 1])

    eQTL <- rbind(eQTL, paste0("cd ", region_dir))
    eQTL <- rbind(eQTL, paste0("cd ", region_dir))
    eQTL <- rbind(eQTL, paste0("mkdir -p ", gene_dir))
    eQTL <- rbind(eQTL, paste0("cd ", gene_dir))

    vcf_dir   <- file.path(INPUT_BASE, paste0("region", region[j, 5]), "temp")
    expr_path <- file.path(region_dir, paste0(region[j, 1], "_expr.txt"))

    eQTL <- rbind(eQTL, paste0(
      ESTEP_MCMC,
      " -vcf ", file.path(vcf_dir, paste0(region[j, 1], ".vcf.gz")),
      " -p ",   expr_path,
      " -maf 0.01 -o eQTL_SumData -LDwindow 1 -GTfield GT -saveSS -zipSS"
    ))
  }
}

write.table(
  eQTL,
  file.path(INPUT_BASE, "eQTL_bash.txt"),
  sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE
)
