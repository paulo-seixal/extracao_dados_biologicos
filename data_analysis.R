library(maftools)
library(DESeq2)
library(ggplot2)
library(summarytools)
library(scales)
library(viridis)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(pheatmap)

gene_data = as.matrix(read.csv('data_stranded_second.csv', row.names = 1))
bio_data = read.csv('bio_data.csv', row.names = 1)
mutation_data = read.csv('mutation_data.csv')


#Seleção de variáveis

definition = bio_data$definition
stage = bio_data$ajcc_clinical_stage
tissue = bio_data$tissue_or_organ_of_origin
alcohol = bio_data$alcohol_history
smoke = bio_data$pack_years_smoked
gender = bio_data$gender
age = bio_data$age_at_index

#Pre processamento -------------------------------------------------------------------------

#tratamento valores omissos
sum(is.na(gene_data))
sum(is.na(mutation_data))
sum(is.na(bio_data))

#retirar amostras de metastases (queremos apenas comparar tecido normal com tecido tumoral primário)
bio_data = bio_data[bio_data$definition != 'Metastatic',]
gene_data = gene_data[,colnames(gene_data) %in% rownames(bio_data)]


#garantir que o barcode das amostras é o mesmo nos dados de gene_data e no bio_data
all(colnames(gene_data) %in% rownames(bio_data)) #False
rownames(bio_data) <- gsub("-", "\\.", rownames(bio_data)) #alterar - por .
all(colnames(gene_data) %in% rownames(bio_data)) # True

#confirmar a ordem dos barcodes das amostras
all(colnames(gene_data) == rownames(bio_data))

#Análise inicial dos dados------------------------------------------------------------------
#tipo de tecido
ggplot(as.data.frame(definition), aes(x = definition, fill = stage)) + 
  geom_bar() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('Frequencias')

#Local do tumor
ggplot(as.data.frame(tissue), aes(x = tissue, fill = stage)) + 
  geom_bar() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('Frequencias')

#Sexo
ggplot(as.data.frame(gender), aes(x = gender)) + 
  geom_bar() + 
  ylab('Frequencias')

#idade
ggplot(as.data.frame(age), aes(x = age)) + 
  geom_histogram() + 
  ylab('Frequencias')

#tabaco

#alcool
#estadiamento

maftools_input = read.maf(mutation_data)

plotmafSummary(maf = maftools_input,
               addStat = 'median',
               dashboard = TRUE)

oncoplot(maf = maftools_input,
         top = 10,
         removeNonMutated = TRUE)

#Analise expressão diferencial--------------------------------------------------------------
dds = DESeqDataSetFromMatrix(countData = gene_data,
                             colData = bio_data,
                             design = ~ definition)


#filter counts under 20
keep = rowSums(counts(dds)) >= 20
dds = dds[keep,]


#run DESeq
dds = DESeq(dds)

#get results
res = results(dds, alpha = 0.05)
diff_exp_analysis <- as.data.frame(res)

#explore results
summary(res)
sum(res$padj < 0.05, na.rm = TRUE)

#filter results
dds_filtered = dds[complete.cases(dds),]


#plotting results

#MA plots
plotMA(res)

#histogram´
hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")



#Enrichment analysis--------------------------------------------------------------------------

#sobre-expressos
genes_up = rownames(res[res$log2FoldChange > 0.5,])
genes_up = gsub("\\..*","",genes_up)

GO_results = enrichGO(gene = genes_up, OrgDb = 'org.Hs.eg.db', keyType = 'ENSEMBL', ont = 'BP', pvalueCutoff = 0.05)

fit = plot(barplot(GO_results, showCategory = 10))


#sub-expressos
genes_down = rownames(res)[res$log2FoldChange < -0.5]
genes_down = gsub("\\..*","",genes_down)

GO_results2 = enrichGO(gene = genes_down, OrgDb = 'org.Hs.eg.db', keyType = 'ENSEMBL', ont = 'BP')

fit2 = plot(barplot(GO_results2, showCategory = 10))


#falta sumarização dos dados
###### estattistica descritiva e gráficos

#falta analise estatistica univariada

#falta analise de expressao diferencial e analise de enriquecimento
