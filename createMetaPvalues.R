# Script takes in the p-values obtained as output of PASTAA, 
# and the correlations of all gene expressions against the 
# gene of interest as inputs and creates a meta-p-value
# resulting in significantly enriched TFs which has
# signficant gene expression correlations

library(metaseqR)

args <- commandArgs(trailingOnly = TRUE)
pastaa_out_file=as.character(args[1])
sc_corr_file=as.character(args[2])
geneOfInterest=toupper(as.character(args[3]))
outfolder=as.character(args[4])

sc_cor <- read.csv(paste0(sc_corr_file,"/",geneOfInterest,"_correlation.csv"), 
                   row.names= 1)
pastaa_out <- read.table(pastaa_out_file, header= F)

cor_pvals <- data.frame(sc_cor_pval=as.numeric(sc_cor[3,]), 
                        row.names = toupper(colnames(sc_cor)))
pastaa_pvals <- data.frame(TF_enrichment_pval= pastaa_out[, 2], 
                            row.names= toupper(pastaa_out[, 1]))
merged_data <- merge(pastaa_pvals, cor_pvals, by= "row.names")
rownames(merged_data) <- merged_data[, 1]
merged_data <- merged_data[, -1]
merged_data[which(merged_data$TF_enrichment_pval > 1), "TF_enrichment_pval"] <- 1

fisher_res <- fisher.method(merged_data)#, na.rm= T)
fisher_res2 <- cbind(fisher_res, merged_data)
fisher_res_sorted <- fisher_res2[order(fisher_res2$p.adj, decreasing= F), ]
fisher_res_sorted$TF_Name<-rownames(fisher_res_sorted)

write.table(fisher_res_sorted, 
            paste0(outfolder,"/Fisher_metapvalues_", geneOfInterest, ".txt"), 
            quote= F, sep= "\t", row.names= F)

