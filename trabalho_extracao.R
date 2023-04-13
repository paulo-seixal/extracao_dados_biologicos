library(TCGAbiolinks)
library(tidyverse)
library(maftools)
library(pheatmap)
library(SummarizedExperiment)
library(GenomicRanges)
library(readr)


#project info
project_summary = getProjectSummary('TCGA-HNSC')



#build query
query = GDCquery(project = 'TCGA-HNSC',
         data.category = 'Transcriptome Profiling',
         experimental.strategy = 'RNA-Seq',
         workflow.type = 'STAR - Counts',
         access = 'open')

query_results = getResults(query)


#download data
GDCdownload(query)


#prepare data
tcga_hnsc = GDCprepare(query, summarizedExperiment = TRUE)
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

#create metadata
bio_data = colData(tcga_hnsc)
bio_data = as.data.frame(bio_data)
write_csv(bio_data, 'bio_data.csv')
