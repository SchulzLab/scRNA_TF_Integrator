# Script takes in a Seurat object and
# calculates correlations of all genes against a gene of interest
# Can perform filtering based on celltype of interest

library(Seurat)

argsexpr <- commandArgs(trailingOnly= T)

defined_args <- c("-file", "-RNA", "-cellannot" ,"-celltype", "-gene", "-output")
arg_tokens <- unlist(strsplit(argsexpr, split= "="))

file_hit <- which(arg_tokens == defined_args[1])
RNA_hit <- which(arg_tokens == defined_args[2])
cellannot_hit <- which(arg_tokens == defined_args[3])
celltype_hit <- which(arg_tokens == defined_args[4])
gene_hit <- which(arg_tokens == defined_args[5])
out_hit <- which(arg_tokens == defined_args[6])

file_name <- arg_tokens[file_hit + 1]
RNA_count_expression <- arg_tokens[RNA_hit + 1]
cellannot <- arg_tokens[cellannot_hit + 1]
celltype <- toupper(arg_tokens[celltype_hit + 1])
geneOfInterest <- toupper(arg_tokens[gene_hit + 1])
outputLoc <- arg_tokens[out_hit + 1]

# read file - can handle only file with Seurat object ().Robj or RDS file)
if(strsplit(file_name, "\\.")[[1]][2] == "rds"){
  print("Reading the Seurat object")
  Sub <- readRDS(file_name)
  Rds_name <- "Sub"
}else{
  print("Loading the Seurat object")
  Rds_name <- load(file_name)
  Sub <- get(Rds_name) ##To make it generalizable for Seurat object stored with any name
}

# Get RNA count data
RNAcounts <- eval(parse(text= paste0(Rds_name, "@", RNA_count_expression))) #Sub@assays$RNA@counts

# Check if celltype info is provided, if so get only RNA counts of that
if(length(celltype_hit)){
        cellannotations <- eval(parse(text= paste0(Rds_name, "@", cellannot)))
        cell_annots <- toupper(as.character(cellannotations))
        cell_annots_idx <- which(cell_annots==celltype)
        RNAcounts <- as.matrix(RNAcounts[,cell_annots_idx])
} else {
        RNAcounts <- as.matrix(RNAcounts)
}
rownames(RNAcounts) <- toupper(rownames(RNAcounts))

# Calculate correlations for the gene of interest against all other genes
print("Calculating correlations of all genes against gene of interest...")
print(dim(RNAcounts))
geneOfInterest_cor <- sapply(seq(nrow(RNAcounts)), function(i) cor.test(RNAcounts[i, ], RNAcounts[geneOfInterest, ]))
colnames(geneOfInterest_cor) <- rownames(RNAcounts)
write.csv(geneOfInterest_cor[seq(7), ], paste0(outputLoc,"/",geneOfInterest,"_correlation.csv"))
print("Saved correlations")
