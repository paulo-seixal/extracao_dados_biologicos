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
library(pheatmap)
library(Glimma)
library(caret)
library(rpart)
library(MLeval)
library(doParallel)
library(RColorBrewer)
library(rattle)

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
head(res.df) #changed ensembl to symbol genes

EnhancedVolcano(res.df, x= 'log2FoldChange', y = 'padj', lab = res.df$symbol,
                pCutoff = 1e-4, FCcutoff = 1) #plot


#MA plots
plotMA(sigs)



#Enrichment analysis--------------------------------------------------------------------------

#sobre-expressos
genes_up = rownames(sigs)[sigs$log2FoldChange > 2]
genes_up = gsub("\\..*","",genes_up)

GO_results = enrichGO(gene = genes_up, OrgDb = 'org.Hs.eg.db', keyType = 'ENSEMBL', ont = 'BP')


#sub-expressos
genes_down = rownames(sigs)[sigs$log2FoldChange < -2]
genes_down = gsub("\\..*","",genes_down)

GO_results2 = enrichGO(gene = genes_down, OrgDb = 'org.Hs.eg.db', keyType = 'ENSEMBL', ont = 'BP')

#plotting
barplot(GO_results, showCategory = 10, title = 'Upregulated')
barplot(GO_results2, showCategory = 10, title = 'Downregulated')

#Filtering dds - apresentar apenas genes significativos (p_adjusted < 0.05)------------------
dds_filtered <- dds[rownames(dds) %in% rownames(sigs), ]
all(rownames(dds_filtered) %in% rownames(sigs)) #True

#Clustering--------------------------------------------------------------
rank = order(sigs$padj)
top_genes <- dds_filtered[rank[1:20], ]
gene_names <- rownames(dds_filtered)  # Extrair nomes dos genes


normal_samples <- top_genes[, colData(top_genes)$definition == "Solid Tissue Normal"]
set.seed(42)  # Set a seed for reproducibility
normal_samples <- normal_samples[, sample(ncol(normal_samples), size = 20)]
tumoral_samples <- top_genes[, colData(top_genes)$definition == "Primary solid Tumor"]
set.seed(42)  # Set a seed for reproducibility
tumoral_samples <- tumoral_samples[, sample(ncol(tumoral_samples), size = 20)]
combined_samples <- cbind(normal_samples, tumoral_samples)

# Calcular distância euclidiana

eucD4 <- dist(t(assay(combined_samples)))
cl.hier5 <- hclust(eucD4)
plot(cl.hier5, labels=colData(combined_samples)$definition)

cl.hier2 <- hclust(eucD4, method="single") 
plot(cl.hier2,labels=colData(combined_samples)$definition)

cl.hier3 <- hclust(eucD4, method="average")
plot(cl.hier3,labels=colData(combined_samples)$definition)
assay_top_genes = t(assay(combined_samples))
ofs = c()

for (k in 2:10) {
  kmeans = kmeans(assay_top_genes, centers = k, nstart = 10)
  ofs = c(ofs, kmeans$tot.withinss)
}

plot = data.frame(num_clusters = 2:10, wss = ofs)

ggplot(plot, aes(x = num_clusters, y = wss)) +
  geom_line() +
  geom_point() +
  labs(x = "Clusters", y = "WSS") +
  ggtitle('Elbow Method')



var_kmeans = kmeans(assay_top_genes, centers = 2)   
var_kmeans
plot(rlog(assay_top_genes), col = var_kmeans$cluster, pch = 16)




#Multidimensional Scaling---------------------------------------------------------------------

glimmaMDS(dds_filtered)


#Machine Learning----------------------------------------------------------------------------

#test and training set
count_data <- t(assay(dds_filtered))
labels <- as.factor(bio_data$definition)

# Split  data
set.seed(123)  # For reproducibility
train_indices <- sample(1:nrow(count_data), nrow(count_data) * 0.75)  # 75% for training
x_train <- count_data[train_indices, ]
y_train <- labels[train_indices]
x_test <- count_data[-train_indices, ]
y_test <- labels[-train_indices]

#transform labels into viable names
y_train = make.names(c(y_train), unique = FALSE)
y_test = make.names(c(y_test), unique = FALSE)

#check sizes
dim(x_train)
length(y_train)
dim(x_test)
length(y_test)

#Função auxiliar para plot da matriz de confusão----------------------------------------


draw_confusion_matrix <- function(cm, title) {
  
  layout(matrix(c(1,1,2)))
  par(mar=c(2,2,2,2))
  plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
  title(title, cex.main=2)
  
  # create the matrix 
  rect(150, 430, 240, 370, col='#3F97D0')
  text(195, 435, 'Solid.Tissue.Normal', cex=1.2)
  rect(250, 430, 340, 370, col='#F7AD50')
  text(295, 435, 'Primary.solid.Tumor', cex=1.2)
  text(125, 370, 'Actual', cex=1.3, srt=90, font=2)
  text(245, 450, 'Predicted', cex=1.3, font=2)
  rect(150, 305, 240, 365, col='#F7AD50')
  rect(250, 305, 340, 365, col='#3F97D0')
  text(140, 400, 'Solid.Tissue.Normal', cex=1.2, srt=90)
  text(140, 335, 'Primary.solid.Tumor', cex=1.2, srt=90)
  
  # add in the cm results 
  res <- as.numeric(cm$table)
  text(195, 400, res[1], cex=1.6, font=2, col='white')
  text(195, 335, res[2], cex=1.6, font=2, col='white')
  text(295, 400, res[3], cex=1.6, font=2, col='white')
  text(295, 335, res[4], cex=1.6, font=2, col='white')
  
  # add in the specifics 
  plot(c(100, 0), c(100, 0), type = "n", xlab="", ylab="", main = "DETAILS", xaxt='n', yaxt='n')
  text(10, 85, names(cm$byClass[1]), cex=1.2, font=2)
  text(10, 70, round(as.numeric(cm$byClass[1]), 3), cex=1.2)
  text(30, 85, names(cm$byClass[2]), cex=1.2, font=2)
  text(30, 70, round(as.numeric(cm$byClass[2]), 3), cex=1.2)
  text(50, 85, names(cm$byClass[5]), cex=1.2, font=2)
  text(50, 70, round(as.numeric(cm$byClass[5]), 3), cex=1.2)
  text(70, 85, names(cm$byClass[6]), cex=1.2, font=2)
  text(70, 70, round(as.numeric(cm$byClass[6]), 3), cex=1.2)
  text(90, 85, names(cm$byClass[7]), cex=1.2, font=2)
  text(90, 70, round(as.numeric(cm$byClass[7]), 3), cex=1.2)
  
  # add in the accuracy information 
  text(30, 35, names(cm$overall[1]), cex=1.5, font=2)
  text(30, 20, round(as.numeric(cm$overall[1]), 3), cex=1.4)
  text(70, 35, names(cm$overall[2]), cex=1.5, font=2)
  text(70, 20, round(as.numeric(cm$overall[2]), 3), cex=1.4)
  
  text(50, -15, "Positive class: Primary.solid.Tumor", cex=1.2, font=2)
} 

#SVM------------------------------------
ctrl = trainControl(method = "repeatedcv", number = 5, repeats = 5, savePredictions = TRUE, classProbs = TRUE)
start.time = proc.time()
svm_model = train(x = x_train, y= y_train, method = 'svmLinear', trControl = ctrl)
stop.time = proc.time()
run.time= stop.time - start.time
print(run.time)

y_svm_predict = predict(svm_model, newdata = x_test)
draw_confusion_matrix(confusionMatrix(as.factor(y_test), y_svm_predict, mode = 'everything'), title = 'Confusion Matrix - SVM')

#SVM com parallel-------------------------------------

cl = makePSOCKcluster(5)
registerDoParallel(cl)
start.time = proc.time()
svm_model = train(x = x_train, y= y_train, method = 'svmLinear', trControl = ctrl)
stop.time = proc.time()
run.time= stop.time - start.time
print(run.time)
stopCluster(cl)

y_svm_predict = predict(svm_model, newdata = x_test)
confusionMatrix(as.factor(y_test), y_svm_predict, mode = 'everything')

importance_svm = varImp(svm_model)
plot(importance_svm, top=20, main='Features (genes) mais importantes - SVM')


#SVM com otimização-----------------------------------------------

# svmGrid <- expand.grid(C=c(0.05, 0.5, 1, 2, 5, 8, 10))
# cl = makePSOCKcluster(5)
# registerDoParallel(cl)
# start.time = proc.time()
# svm_model_tuned = train(x = x_train, y= y_train, method = 'svmLinear', trControl = ctrl, tuneGrid = svmGrid)
# stop.time = proc.time()
# run.time= stop.time - start.time
# print(run.time)
# stopCluster(cl)
# y_svm_grid_predict = predict(svm_model_tuned, newdata = x_test)
# confusionMatrix(as.factor(y_test), y_svm_grid_predict, mode = 'everything')

#Decision Tree---------------------------------

cl = makePSOCKcluster(5)
registerDoParallel(cl)
start.time = proc.time()
decision_tree_model = train(x = x_train, y= y_train, method = 'rpart', trControl = ctrl)
stop.time = proc.time()
run.time= stop.time - start.time
print(run.time)
stopCluster(cl)
y_tree_predict = predict(decision_tree_model, newdata = x_test)
confusionMatrix(as.factor(y_test), y_tree_predict, mode = 'everything')

importance_tree= varImp(decision_tree_model)
plot(importance_tree, top=20, main='Features (genes) mais importantes - Decision Tree')

fancyRpartPlot(decision_tree_model$finalModel)
#MLP------------------------------------

cl = makePSOCKcluster(5)
registerDoParallel(cl)
start.time = proc.time()
mlp_model = train(x = x_train, y= y_train, method = 'mlp', trControl = ctrl)
stop.time = proc.time()
run.time= stop.time - start.time
print(run.time)
stopCluster(cl)
y_mlp_predict = predict(mlp_model, newdata = x_test)
draw_confusion_matrix(confusionMatrix(as.factor(y_test), y_mlp_predict, mode = 'everything'), title = 'Confusion Matrix - MLP')

importance_mlp = varImp(mlp_model)
plot(importance_mlp, top=20, main='Features (genes) mais importantes - Multi-Layer Perceptron')



#KNN---------------------------------------------------
cl = makePSOCKcluster(5)
registerDoParallel(cl)
start.time = proc.time()
knn_model = train(x = x_train, y= y_train, method = 'knn', trControl = ctrl)
stop.time = proc.time()
run.time= stop.time - start.time
print(run.time)
stopCluster(cl)
y_knn_predict = predict(knn_model, newdata = x_test)
confusionMatrix(as.factor(y_test), y_knn_predict, mode = 'everything')

importance_knn = varImp(knn_model)
plot(importance_knn, top=20, main='Features (genes) mais importantes - KNN')
#KNN com otimização--------------------------------------
#knnGrid = expand.grid(k =c(1:10) )
#cl = makePSOCKcluster(5)
#registerDoParallel(cl)
#start.time = proc.time()
#knn_model_tuned = train(x = x_train, y= y_train, method = 'knn', trControl = ctrl, tuneGrid = knnGrid)
#stop.time = proc.time()
#run.time= stop.time - start.time
#print(run.time)
#stopCluster(cl)
#y_knn_grid_predict = predict(knn_model_tuned, newdata = x_test)
#confusionMatrix(as.factor(y_test), y_knn_grid_predict, mode = 'everything')

#Naive Bayes---------------------------------------------------
cl = makePSOCKcluster(5)
registerDoParallel(cl)
start.time = proc.time()
nb_model = train(x = x_train, y= y_train, method = 'naive_bayes', trControl = ctrl)
stop.time = proc.time()
run.time= stop.time - start.time
print(run.time)
stopCluster(cl)
y_nb_predict = predict(nb_model, newdata = x_test)
confusionMatrix(as.factor(y_test), y_nb_predict, mode = 'everything')

importance_nb = varImp(nb_model)
plot(importance_nb, top=20, main='Features (genes) mais importantes - Naive Bayes')

#Curvas ROC------------------------
res_roc = evalm(list(decision_tree_model, knn_model, mlp_model, nb_model, svm_model), positive = 'Primary.solid.Tumor', gnames = c('Decision Tree','KNN','Multilayer Perceptron', 'Naive Bayes','SVM'), plots = 'r', title = 'ROC-Curves', rlinethick = 0.8, cols = c('red','blue','green','darkgreen','purple'))
res_pr = evalm(list(decision_tree_model, knn_model, mlp_model,nb_model, svm_model), positive = 'Primary.solid.Tumor', gnames = c('Decision Tree','KNN','Multilayer Perceptron','Naive Bayes', 'SVM'), plots = 'prg', title = 'Precision-Recall Curves', rlinethick = 0.8, cols = c('red','blue','green','darkgreen','purple'))

y_svm_predict_prob = predict(svm_model, newdata = x_test, type = 'prob')
y_tree_predict_prob = predict(decision_tree_model, newdata = x_test, type = 'prob')
y_mlp_predict_prob = predict(mlp_model, newdata = x_test, type = 'prob')
y_knn_predict_prob = predict(knn_model, newdata = x_test, type = 'prob')
y_nb_predict_prob = predict(nb_model, newdata = x_test, type = 'prob')

y_svm_predict_prob$obs = y_test
y_svm_predict_prob$Group = 'SVM'
y_tree_predict_prob$obs = y_test
y_tree_predict_prob$Group = 'Decision Tree'
y_mlp_predict_prob$obs = y_test
y_mlp_predict_prob$Group = 'Multilayer Perceptron'
y_knn_predict_prob$obs = y_test
y_knn_predict_prob$Group = 'KNN'
y_nb_predict_prob$obs = y_test
y_nb_predict_prob$Group = 'Naive Bayes'

combo_df = rbind(y_svm_predict_prob, y_tree_predict_prob, y_mlp_predict_prob, y_knn_predict_prob, y_nb_predict_prob)

test_roc = evalm(combo_df, plots = 'r', title = 'ROC-Curves', rlinethick = 0.8, cols = c('red','blue','green','darkgreen','purple'))
test_pr = evalm(combo_df, plots = 'prg', title = 'Precision-Recall Curves', rlinethick = 0.8, cols = c('red','blue','green','darkgreen','purple'))
