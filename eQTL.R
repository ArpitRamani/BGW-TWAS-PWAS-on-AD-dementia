options(scipen = 999)
expr <- read.delim("/home/jyang51/YangLabData/jyang/BGW_TPWAS_AD/Data/RNAseq_BulkBrain/ROSMAP_expr_TIGAR_format_WGS_IDs_b38_2023_chr1-22.tsv")
regions <- read.delim(paste0("/home/jyang51/YangLabData/aramani/GIFT2/input1/region_from_matched.txt",sep=""))
for (i in 1:length(unique(regions[,"region"]))){
  region <- regions[which(regions[,"region"]==i),]
  region <- region[order(region$GeneName), ]
  for (j in 1:nrow(region)){
    gene_expr <- expr[which(expr[,"GeneName"]==region[j,"GeneName"]),]
    gene_expr <- t(gene_expr[,6:ncol(gene_expr)])
    gene_expr <- cbind(rownames(gene_expr),gene_expr)
    gene_expr[,1] <- gsub("\\.","-",gene_expr[,1])
    gene_expr <- na.omit(gene_expr)
    write.table(gene_expr,paste0("/home/jyang51/YangLabData/aramani/GIFT2/input1/region",
                                   i,"/",region[j,1],"_expr.txt",sep=""),sep="\t",quote = FALSE,col.names = FALSE,row.names = FALSE)
  }
}

options(scipen = 999)
eQTL <- c()
regions <- read.delim(paste0("/home/jyang51/YangLabData/aramani/GIFT2/input1/region_from_matched.txt",sep=""))
for (i in 1:length(unique(regions[,"region"]))){
    region <- regions[which(regions[,"region"]==i),]
    region <- region[order(region$GeneName), ]
    for (j in 1:nrow(region)){
      eQTL <- rbind(eQTL,paste0("cd /home/jyang51/YangLabData/aramani/GIFT2/input1/region",i,sep=""))
      eQTL <- rbind(eQTL,paste0("cd /home/jyang51/YangLabData/aramani/GIFT2/input1/region",i,sep=""))
      eQTL <- rbind(eQTL,paste0("mkdir -p /home/jyang51/YangLabData/aramani/GIFT2/input1/region",i,"/",region[j,1],sep=""))
      eQTL <- rbind(eQTL,paste0("cd /home/jyang51/YangLabData/aramani/GIFT2/input1/region",i,"/",region[j,1],sep=""))

      eQTL <- rbind(eQTL,paste0(
        "/home/jyang51/jyang/GITHUB/BGW-xWAS-SS/bin/Estep_mcmc -vcf /home/jyang51/YangLabData/aramani/GIFT2/input1/region",
        region[j,5], "/temp/", region[j,1], ".vcf.gz -p /home/jyang51/YangLabData/aramani/GIFT2/input1/region",
        i, "/", region[j,1], "_expr.txt -maf 0.01 -o eQTL_SumData -LDwindow 1 -GTfield GT -saveSS -zipSS",
        sep=""
      ))
  }
}
write.table(eQTL,paste0("/home/jyang51/YangLabData/aramani/GIFT2/input1/eQTL_bash.txt",sep=""),
            sep="\t",quote = FALSE,col.names = FALSE,row.names = FALSE)
