# Script takes in MetaCellaR created metacell based scRNA RNAcounts
# and calculates correlations of all genes against a gene of interest

argsexpr <- commandArgs(trailingOnly= T)

defined_args <- c("-file" ,"-celltype", "-gene", "-output")
arg_tokens <- unlist(strsplit(argsexpr, split= "="))

file_hit <- which(arg_tokens == defined_args[1])
celltype_hit <- which(arg_tokens == defined_args[2])
gene_hit <- which(arg_tokens == defined_args[3])
out_hit <- which(arg_tokens == defined_args[4])

file_name <- arg_tokens[file_hit + 1]
geneOfInterest <- toupper(arg_tokens[gene_hit + 1])
outputLoc <- arg_tokens[out_hit + 1]

df=read.csv(file_name, header = T, stringsAsFactors = F)
#df=read.csv("./cellSummarized_kmed_means.csv", header = T, stringsAsFactors = F)
#cellOfInterest="kidney proximal straight tubule epithelial cell"
#geneOfInterest <- toupper("Nox4")
if(length(celltype_hit)){
  cellOfInterest <- arg_tokens[celltype_hit + 1]
  # replace spaces in the celltype name
  cellOfInterestMod=gsub(" ",".", cellOfInterest)
  cellidx=grep(cellOfInterestMod, colnames(df), ignore.case = T)
  if(length(cellidx)>0){
    df_cells=df[,c(1,cellidx)]
    print("subsetted the celltype of interest...")
  } else {
    df_cells=df
    print("Unable to find celltype of interest, using all celltypes available...")
  }
} else {
  df_cells=df
  print("No celltype of interest provided, using all celltypes available...")
}
print(c("Before celltype subset: ", dim(df)))
print(c("After celltype subset:", dim(df_cells)))

rownames(df_cells)=toupper(df_cells[,1])
df_cells_m=as.matrix(df_cells[,-1])

# Calculate correlations for the gene of interest against all other genes
print("Calculating correlations of all genes against gene of interest...")
geneOfInterest_cor <- sapply(seq(nrow(df_cells)),
                             function(i) cor.test(df_cells_m[i, ],
                                                  df_cells_m[geneOfInterest, ]))
colnames(geneOfInterest_cor) <- rownames(df_cells)
write.csv(geneOfInterest_cor[seq(7), ], paste0(outputLoc,"/",geneOfInterest,"_correlation.csv"))
print("Saved correlations")
