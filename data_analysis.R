library(maftools)
library(DESeq2)
gene_data = read.csv('data_stranded_second.csv', row.names = 1)
bio_data = read.csv('bio_data.csv', row.names = 1)
mutation_data = read.csv('mutation_data.csv')


#Pre processamento -------------------------------------------------------------------------

#tratamento valores omissos
sum(is.na(gene_data))
sum(is.na(mutation_data))
sum(is.na(bio_data))

#filtrar contagens abixo de 10
keep = rowSums(counts(gene_data)) >= 10
gene_data = gene_data[keep,]

#garantir que o barcode das amostras é o mesmo nos dados de gene_data e no bio_data
all(colnames(gene_data) %in% rownames(bio_data)) #False
colnames(gene_data) <- gsub("\\.", "-", colnames(gene_data)) #alterar . por -
all(colnames(gene_data) %in% rownames(bio_data)) # True

#confirmar a ordem dos barcodes das amostras
all(colnames(gene_data) == rownames(bio_data))

#




maftools_input = read.maf(mutation_data)

plotmafSummary(maf = maftools_input,
               addStat = 'median',
               dashboard = TRUE)

oncoplot(maf = maftools_input,
         top = 10,
         removeNonMutated = TRUE)


#falta sumarização dos dados
###### estattistica descritiva e gráficos

#falta analise estatistica univariada

#falta analise de expressao diferencial e analise de enriquecimento
