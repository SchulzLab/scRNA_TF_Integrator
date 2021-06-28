#!/usr/bin/bash

GENE_NAME_OF_INTEREST=""
EPIREGIO_FILE=""
OUTPUT_FOLDER=""
GENOME_SEQ_FASTA="/home/skarunanithi/identifyEnrichedTFsWithscRNAdata/scripts/hg38.fa"
PATH_TO_TF_BINDING_WORKFLOW="/home/skarunanithi/identifyEnrichedTFsWithscRNAdata/scripts"
PATH_TO_METACELLAR="/home/skarunanithi/identifyEnrichedTFsWithscRNAdata/MetaCellaR"
PATH_TO_scRNA_PIPELINE="/home/skarunanithi/identifyEnrichedTFsWithscRNAdata/scripts"
SINGLE_CELL_EXPRESSION_FILE=""
RNA_EXPRESSION_SLOT=""
CELL_TYPE_ANNOTATION_SLOT=""
CELL_TYPE_OF_INTEREST="none"
RUN_UMAP="T"
GENE_NAME_OF_INTEREST=""

function usage()
{
	echo "Usage: "
    echo "bash ~/identifyEnrichedTFsWithscRNAdata/scripts/Pipeline.sh -o=./testNOX4 -e=NOX4_REMs.csv -mcr=/home/skarunanithi/identifyEnrichedTFsWithscRNAdata/MetaCellaR -scf=./droplet_Liver_seurat_tiss.Robj -r='data' -c='meta.data$cell_ontology_class' -u=T -g=NOX4 -cell='hepatocyte'"
    echo "pipeline.sh "
    echo "Parameters:"
    echo "-h|--help"
    echo "-o | --outloc - path to the output directory (will be created if not available)"
    echo "-e | --epiregio - Epiregio ouptut for the gene of your interest downloaded as a CSV file "
    echo "-gseq | --genome_file - genome sequence file in FASTA format (Default: /home/skarunanithi/identifyEnrichedTFsWithscRNAdata/scripts/hg38.fa)"
    echo "-tfw | --tfbinding_workflow_loc - Path to the TRAP/FIMO workflow created by Nina (Default: /home/skarunanithi/identifyEnrichedTFsWithscRNAdata/scripts)"
    echo "-tfw_scr | --tfbinding_scRNA_workflow_loc - Path to scRNA integration pipeline created by Siva, which has additional processing scripts (Default: /home/skarunanithi/identifyEnrichedTFsWithscRNAdata/scripts) "
    echo "-mcr | --metacellar_loc - Path to MetaCellaR created by Fatemeh (Default: /home/skarunanithi/identifyEnrichedTFsWithscRNAdata/MetaCellaR)"
    echo "-scf | --scRNA_file - Name of the scRNA file containing a Seurat object"
    echo "-r | --scRNA_slot - Name of the RNA expression data slot in the Seurat object"
    echo "-c | --scRNA_cell_annot_slot - Name of the cell type annotation data slot in the Seurat object"
    echo "-u | --run_umap - Boolean T or F, whether to run UMAP as part of MetaCellaR (Default: T)"
    echo "-g | --gene_name - Name of the gene of interest (Eg. NOX4)"
    echo "-cell | --cell_name - Name of the interested cell type (if any; Default: None)"
}

###### parameter parsing ######
while [ "$1" != "" ]; do
    PARAM=`echo $1 | awk -F= '{print $1}'`
    VALUE=`echo $1 | awk -F= '{print $2}'`
    case $PARAM in
        -h | --help)
            usage
            exit
            ;;
        -o | --outloc)
            OUTPUT_FOLDER=$VALUE
            ;;
        -e | --epiregio)
            EPIREGIO_FILE=$VALUE
            ;;
        -gseq | --genome_file)
            GENOME_SEQ_FASTA=$VALUE
            ;;
        -tfw | --tfbinding_workflow_loc)
            PATH_TO_TF_BINDING_WORKFLOW=$VALUE
            ;;
        -tfw_scr | --tfbinding_scRNA_workflow_loc)
            PATH_TO_scRNA_PIPELINE=$VALUE
            ;;
        -mcr | --metacellar_loc)
            PATH_TO_METACELLAR=$VALUE
            ;;
        -scf | --scRNA_file)
            SINGLE_CELL_EXPRESSION_FILE=$VALUE
            ;;
        -r | --scRNA_slot)
            RNA_EXPRESSION_SLOT=$VALUE
            ;;
        -c | --scRNA_cell_annot_slot)
            CELL_TYPE_ANNOTATION_SLOT=$VALUE
            ;;
        -u | --run_umap)
            RUN_UMAP=$VALUE
            ;;
        -g | --gene_name)
            GENE_NAME_OF_INTEREST=$VALUE
            ;;
        -cell | --cell_name)
            CELL_TYPE_OF_INTEREST=$VALUE
            ;;
        *)
            echo "ERROR: unknown parameter \"$PARAM\""
            usage
            exit 1
            ;;
    esac
    shift
done
###### Done parameter parsing ######

METACELLAR_OUTPUT=${OUTPUT_FOLDER}/metacellar_out
PASTAA_OUTPUT=${OUTPUT_FOLDER}/pastaa_out
FIMO_OUTPUT=${OUTPUT_FOLDER}/fimo_out

###### TO_DO: Error handling #######

###### Nina's scripts to identify enriched TF Binding Sites using PASTAA and FIMO  ######

bash ${PATH_TO_TF_BINDING_WORKFLOW}/PASTAA_Fimo_analysis/workflow.sh ${EPIREGIO_FILE} ${FIMO_OUTPUT}/ ${PASTAA_OUTPUT}/ ${GENOME_SEQ_FASTA} ${PATH_TO_TF_BINDING_WORKFLOW}/

echo "TF Binding Site Prediction Workflow Completed"

###### Fatemeh's MetaCellaR run #######

#if [ ${CELL_TYPE_ANNOTATION_SLOT} == "" ]
#then
#  Rscript ${PATH_TO_METACELLAR}/MetaCellaR.R -file ${SINGLE_CELL_EXPRESSION_FILE} -RNA ${RNA_EXPRESSION_SLOT} -umap ${RUN_UMAP} -output ${METACELLAR_OUTPUT}
#else
Rscript ${PATH_TO_METACELLAR}/MetaCellaR.R -file ${SINGLE_CELL_EXPRESSION_FILE} -RNA ${RNA_EXPRESSION_SLOT} -celltype ${CELL_TYPE_ANNOTATION_SLOT} -umap ${RUN_UMAP} -output ${METACELLAR_OUTPUT}

echo "Identification of Metacells Completed"

###### Calculate correlation of single cell (metacell) gene expression for all genes against gene of interest ######

if [ $CELL_TYPE_OF_INTEREST == "none" ]; then
  Rscript ${PATH_TO_scRNA_PIPELINE}/createMetaCellGeneExpCorrelations.R -file=${METACELLAR_OUTPUT}/cellSummarized_kmed_means.csv -gene=${GENE_NAME_OF_INTEREST} -output=${OUTPUT_FOLDER}
else
  Rscript ${PATH_TO_scRNA_PIPELINE}/createMetaCellGeneExpCorrelations.R -file=${METACELLAR_OUTPUT}/cellSummarized_kmed_means.csv -celltype="$CELL_TYPE_OF_INTEREST" -gene=${GENE_NAME_OF_INTEREST} -output=${OUTPUT_FOLDER}
fi

###### Creating meta p-value from the outputs of MetaCellaR and TF-Enrichment ######
echo "Integrating single cell and TF enrichment results... "
Rscript ${PATH_TO_scRNA_PIPELINE}/createMetaPvalues.R ${PASTAA_OUTPUT}/PASTAA_output.txt ${OUTPUT_FOLDER} ${GENE_NAME_OF_INTEREST} ${OUTPUT_FOLDER}

echo "Integrating single cell and TF enrichment results Completed. Results can be found in "${OUTPUT_FOLDER}
