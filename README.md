# scRNA_TF_Integrator
Integrate scRNA data with TF binding site prediction analysis. These scripts are still under development and are currently being used only internally by the [SchulzLab](https://schulzlab.github.io/) group members. 

The workflow is as follows.

 - Step-1: For a gene of interest get the REMs from the EpiRegioDB as a CSV file
 - Step-2: We run the scripts to identify enriched TFs in the TF Binding Sites which can be found in the REMs using PASTAA and FIMO (created by Nina). This step gives a list of enriched TFs and their p-values. The files associated with this step are located in scripts/src, scripts/PASTAA_Fimo_analysis, scripts/identifyEnrichedTFs, and scripts/identifyTFBindingSites.
 - Step-3: We run MetaCellaR (created by Fatemeh) to create Metacells from the available single cell data. The current version of MetaCellaR can be found at https://github.com/SchulzLab/MetaCellaR.
 - Step-4: Using the Metacells identified in Step-3, we calculate the RNA expression correlation of your gene of interest against all genes, and their p-values. You can subset this step to include only a particular cell type of your interest. If the interested cell type is unavailable, all available cell types will be used.
 - Step-5: We calculate a meta p-value (uses the metaseqR package) by integrating the p-values from Step-2 and Step-4. Please note that this step will report the meta p-values for only the enriched TFs identified in Step-2.

## Example run

```bash ~/identifyEnrichedTFsWithscRNAdata/scripts/Pipeline.sh -o=./testNOX4 -e=NOX4_REMs.csv -mcr=~/identifyEnrichedTFsWithscRNAdata/MetaCellaR -scf=./droplet_Liver_seurat_tiss.Robj -r='data' -c='meta.data$cell_ontology_class' -u=T -g=NOX4 -cell='hepatocyte'```

## Description of commands


| Short        | Long           | Description  |
|:-------------|:-------------|:-----|
| -h | --help | Display Help commands |
| -o | --outloc |Path to the output directory (will be created if not available|
|-e|--epiregio | Epiregio ouptut for the gene of your interest downloaded as a CSV file |
|-gseq | --genome_file | genome sequence file in FASTA format (Default: ~/identifyEnrichedTFsWithscRNAdata/scripts/hg38.fa) |
|-tfw | --tfbinding_workflow_loc | Path to the TRAP/FIMO workflow created by Nina (Default: ~/identifyEnrichedTFsWithscRNAdata/scripts)|
|-tfw_scr | --tfbinding_scRNA_workflow_loc | Path to scRNA integration pipeline created by Siva, which has additional processing scripts (Default: ~/identifyEnrichedTFsWithscRNAdata/scripts)|
|-mcr | --metacellar_loc | Path to MetaCellaR created by Fatemeh (Default: ~/identifyEnrichedTFsWithscRNAdata/MetaCellaR)|
|-scf | --scRNA_file | Name of the scRNA file containing a Seurat object|
|-r | --scRNA_slot | Name of the RNA expression data slot in the Seurat object|
|-c | --scRNA_cell_annot_slot | Name of the cell type annotation data slot in the Seurat object|
|-u | --run_umap | Boolean T or F, whether to run UMAP as part of MetaCellaR (Default: T)|
|-g | --gene_name | Name of the gene of interest (Eg. NOX4)|
|-cell | --cell_name | Name of the interested cell type (if any; Default: None)|
