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
library(haven)
library(textshaping)
library(ggrepel)
library(EnhancedVolcano)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(genefilter)


#import
gene_data = as.matrix(read.csv('data_unstranded.csv', row.names = 1))
gene_data_fpkm = as.matrix(read.csv('data_fpkm_unstrand.csv', row.names = 1))
bio_data = read.csv('bio_data.csv', row.names = 1)
mutation_data = read.csv('mutation_data.csv')

#Seleção de variáveis --------------------------------------------------------------------

definition = bio_data$definition
stage = bio_data$ajcc_clinical_stage
tissue = bio_data$tissue_or_organ_of_origin
alcohol = bio_data$alcohol_history
smoke = bio_data$pack_years_smoked
smoke_disc = ifelse(is.na(smoke), "No", ifelse(smoke > 0, "Yes", "No"))
sum(is.na(smoke)) == sum(smoke_disc=='No')
gender = bio_data$gender
age = bio_data$age_at_index

#Pre processamento -------------------------------------------------------------------------

#tratamento valores omissos
sum(is.na(gene_data))
sum(is.na(mutation_data))
sum(is.na(bio_data))

#garantir que o barcode das amostras é o mesmo nos dados de gene_data e no bio_data
all(colnames(gene_data) %in% rownames(bio_data)) #False
rownames(bio_data) <- gsub("-", "\\.", rownames(bio_data)) #alterar - por .
all(colnames(gene_data) %in% rownames(bio_data)) # True

#retirar amostras de metastases (queremos apenas comparar tecido normal com tecido tumoral primário no DESeq2)
bio_data = bio_data[bio_data$definition != 'Metastatic',]
gene_data = gene_data[,colnames(gene_data) %in% rownames(bio_data)]
gene_data_fpkm = gene_data_fpkm[,colnames(gene_data_fpkm) %in% rownames(bio_data)]
all(colnames(gene_data) %in% rownames(bio_data)) # True

#confirmar a ordem dos barcodes das amostras
all(colnames(gene_data) == rownames(bio_data))

#realizar t-test e verificar p-values e reirar os genes com maior evidência de expressão diferencial
tt = rowttests(gene_data_fpkm, as.factor(bio_data$definition))
rank = order(tt$p.value)
p10000 = rank[1:10000]
rows = which(tt$p.value %in% tt$p.value[p10000])
rownames(tt)[rows] #genes com menor p-value, logo com maior probabilidade de serem diferencialmente expressos

#filtrar dataset de contagens v«elos valores anteriores
gene_data = gene_data[rownames(gene_data) %in% rownames(tt)[rows],]
dim(gene_data) #passa a contar apenas com os 10000 genes de menor p-value


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
ggplot(as.data.frame(gender), aes(x = gender, fill = gender)) + 
  geom_bar() +
  geom_text(aes(label = paste0(round((..count..)/sum(..count..) * 100), "%")), 
            stat = "count", vjust = -0.5, size = 5) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('Frequencias')

#idade
ggplot(as.data.frame(age), aes(x = age)) + 
  geom_histogram(bins = 20, fill = 'lightblue', color = 'black') + 
  labs(title = 'Idades', x = 'Idade', y = 'Frequencias')

#tabaco
ggplot(as.data.frame(smoke_disc), aes(x = smoke_disc)) +
  geom_bar() +
  ylab('Frequencias')

#alcool
ggplot(as.data.frame(alcohol), aes(x = alcohol, fill = alcohol)) + 
  geom_bar() +
  geom_text(aes(label = paste0(round((..count..)/sum(..count..) * 100), "%")), 
            stat = "count", vjust = -0.5, size = 5) + 
  ylab('Frequencies')

#estadiamento
ggplot(as.data.frame(stage),aes(x = stage, fill= stage)) +
  geom_bar() +
  ylab('Frequencias')

s = as.data.frame(stage)
stage_freq = table(s$stage)
df = data.frame(stage = stage(stage_freq), freq = as.numeric(stage_freq))

slices = df$freq
names = df$stage.Var1
pct = round(slices/sum(slices)*100)

new_labels = paste(names, ' - ', slices, '(', pct, '%',')', sep="")

pie(slices, labels = new_labels, main = "Estadiamento do tumor", col = rainbow(6))

#Análise Sumária dos dados ------------------------------------------------------------------------
view(dfSummary(bio_data))
view(dfSummary(gene_data))
view(dfSummary(mutation_data))

#Mutações ----------------------------------------------------------------------------------------

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


#set reference level
bio_data$definition = relevel(as.factor(bio_data$definition), ref = 'Solid Tissue Normal')


#run DESeq
dds = DESeq(dds)


#Quality control - idealmente condições semelhantes no mesmo cluster
vsdata = vst(dds, blind = FALSE)
plotPCA(vsdata, intgroup = 'definition')

#dispersion plot
plotDispEsts(dds)

#get results from DESeq 
res = results(dds, alpha = 0.05)
summary(res)

#signifiatives - only padj under 0.05
sigs = na.omit(res)
sigs = sigs[sigs$padj < 0.05,]

#explore results
summary(sigs)

#volcano plot
res.df = as.data.frame(res)
rownames(res.df) = gsub("\\..*","",rownames(res.df))
res.df$symbol = mapIds(org.Hs.eg.db, keys = rownames(res.df), keytype = 'ENSEMBL', column = 'SYMBOL')
res.df #changed ensembl to symbol genes

EnhancedVolcano(res.df, x= 'log2FoldChange', y = 'padj', lab = res.df$symbol,
                pCutoff = 1e-4, FCcutoff = 1) #plot


#MA plots
plotMA(sigs)


#heatmap
sigs.df = as.data.frame(sigs)
rownames(sigs.df) = gsub("\\..*","",rownames(sigs.df))
sigs.df$symbol = mapIds(org.Hs.eg.db, keys = rownames(sigs.df), keytype = 'ENSEMBL', column = 'SYMBOL')
sigs.df #changed ensembl to symbol genes

sigs.df_top = sigs.df[(sigs.df$baseMean > 50) & (abs(sigs.df$log2FoldChange) > 5), ] #filtering for top genes
sigs.df_top = sigs.df_top[order(sigs.df_top$log2FoldChange, decreasing = TRUE),] #order by fold change
sigs.df_top #up regulated first, down regulated last

rlog_out = vst(dds, blind = FALSE) #normalized count from dds
rownames(rlog_out) = gsub("\\..*","",rownames(rlog_out))
mat = assay(rlog_out)[rownames(sigs.df_top), rownames(bio_data)] #matrix sig genes x samples
colnames(mat) = NULL
base_mean = rowMeans(mat)
mat.scaled = t(apply(mat, 1, scale)) #center and scale each column, then transpose
colnames(mat.scaled) = colnames(mat)

num_keep = 25
rows_keep = c(seq(1:num_keep), seq((nrow(mat.scaled)-num_keep), nrow(mat.scaled)))

l2_val =as.matrix(sigs.df_top[rows_keep,]$log2FoldChange) 
colnames(l2_val) = 'logFC'

col_logFC = colorRamp2(c(min(l2_val),0, max(l2_val)), c("blue", "white", "red")) 



ha <- HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = 2), 
                                               height = unit(2, "cm")))
h1 <- Heatmap(mat.scaled[rows_keep,], cluster_rows = F, 
              column_labels = colnames(mat.scaled), name="Z-score",
              cluster_columns = T)
h2 <- Heatmap(l2_val, row_labels = sigs.df_top$symbol[rows_keep], 
              cluster_rows = F, name="logFC", top_annotation = ha, col = col_logFC)

h<-h1+h2
h


#Enrichment analysis--------------------------------------------------------------------------

#sobre-expressos
genes_up = rownames(sigs.df)[sigs.df$log2FoldChange > 2]
genes_up = gsub("\\..*","",genes_up)

GO_results = enrichGO(gene = genes_up, OrgDb = 'org.Hs.eg.db', keyType = 'ENSEMBL', ont = 'BP')


#sub-expressos
genes_down = rownames(sigs.df)[sigs.df$log2FoldChange < -2]
genes_down = gsub("\\..*","",genes_down)

GO_results2 = enrichGO(gene = genes_down, OrgDb = 'org.Hs.eg.db', keyType = 'ENSEMBL', ont = 'BP')

#plotting
barplot(GO_results, showCategory = 10, title = 'Upregulated')
barplot(GO_results2, showCategory = 10, title = 'Downregulated')
goplot(GO_results2)

#falta sumarização dos dados
###### estattistica descritiva e gráficos

#falta analise estatistica univariada

#falta analise de expressao diferencial e analise de enriquecimento

