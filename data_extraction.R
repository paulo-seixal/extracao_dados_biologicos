library(TCGAbiolinks)
library(tidyverse)
library(maftools)
library(pheatmap)
library(SummarizedExperiment)



#project info
project_summary = getProjectSummary('TCGA-HNSC')



#build query gene expression
query_gene = GDCquery(project = 'TCGA-HNSC',
         data.category = 'Transcriptome Profiling',
         experimental.strategy = 'RNA-Seq',
         workflow.type = 'STAR - Counts',
         access = 'open')

query_gene_results = getResults(query_gene)


#download data
GDCdownload(query_gene)


#prepare data
tcga_hnsc = GDCprepare(query_gene, summarizedExperiment = TRUE)
data_unstranded = assay(tcga_hnsc, 'unstranded')
data_stranded_first = assay(tcga_hnsc, 'stranded_first')
data_stranded_second = assay(tcga_hnsc, 'stranded_second')
data_tpm_unstrand = assay(tcga_hnsc, 'tpm_unstrand')
data_fpkm_unstrand = assay(tcga_hnsc, 'fpkm_unstrand')
data_fpkm_uq_unstrand = assay(tcga_hnsc, 'fpkm_uq_unstrand')

#create files
write.csv(data_unstranded, 'data_unstranded.csv')
write.csv(data_stranded_first, 'data_stranded_first.csv')
write.csv(data_stranded_second, 'data_stranded_second.csv')
write.csv(data_tpm_unstrand, 'data_tpm_unstrand.csv')
write.csv(data_fpkm_unstrand, 'data_fpkm_unstrand.csv')
write.csv(data_fpkm_uq_unstrand, 'data_fpkm_uq_unstrand.csv')

#create biological data
bio_data = colData(tcga_hnsc)
bio_data = as.data.frame(bio_data)
write_csv(bio_data, 'bio_data.csv')


#-------------------------------------------------------------------------------------------#


#build query mutation data
query_mutation = GDCquery(project = 'TCGA-HNSC',
                          data.category = 'Simple Nucleotide Variation',
                          access = 'open')

query_mutation_results = getResults(query_mutation)


#download data
GDCdownload(query_mutation)


#prepare mutation data
mutation_data = GDCprepare(query_mutation, summarizedExperiment = TRUE)

#create files
write_csv(mutation_data, 'mutation_data.csv')

