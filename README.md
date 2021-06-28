# scRNA_TF_Integrator
Integrate scRNA data with TF binding site prediction analysis. These scripts are still under development and are currently being used only internally by the (SchulzLab)[https://schulzlab.github.io/] group members. 

The workflow is as follows.

 - Step-1: For a gene of interest get the REMs from the EpiRegioDB as a CSV file
 - Step-2: We run the scripts to identify enriched TFs in the TF Binding Sites which can be found in the REMs using PASTAA and FIMO (created by Nina). This step gives a list of enriched TFs and their p-values. The files associated with this step are located in scripts/src, scripts/PASTAA_Fimo_analysis, scripts/identifyEnrichedTFs, and scripts/identifyTFBindingSites.
 - Step-3: We run MetaCellaR (created by Fatemeh) to create Metacells from the available single cell data. The MetaCellaR scripts should be available locally in the MetaCellaR folder of this directory. The current version of MetaCellaR can be found at https://github.com/SchulzLab/MetaCellaR.
 - Step-4: Using the Metacells identified in Step-3, we calculate the RNA expression correlation of your gene of interest against all genes, and their p-values. You can subset this step to include only a particular cell type of your interest. If the interested cell type is unavailable, all available cell types will be used.
 - Step-5: We calculate a meta p-value by integrating the p-values from Step-2 and Step-4. Please note that this step will report the meta p-values for only the enriched TFs identified in Step-2.

Example run:

```bash ~/identifyEnrichedTFsWithscRNAdata/scripts/Pipeline.sh -o=./testNOX4 -e=NOX4_REMs.csv -mcr=/home/skarunanithi/identifyEnrichedTFsWithscRNAdata/MetaCellaR -scf=./droplet_Liver_seurat_tiss.Robj -r='data' -c='meta.data$cell_ontology_class' -u=T -g=NOX4 -cell='hepatocyte```

