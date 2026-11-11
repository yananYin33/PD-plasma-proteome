
rm(list=ls())
dev.off()
library(WGCNA)
library(reshape2)
library(stringr)
exprs_df <- read.csv('1.csv',row.names = 1)
exprs_df [1:3,1:4]
dim(exprs_df)
names(exprs_df)
datExpr <- log2(exprs_df +1)
dataExpr <- as.data.frame(t(datExpr))
head(dataExpr)
dim(dataExpr)
table(is.na(dataExpr))
gsg = goodSamplesGenes(dataExpr, verbose = 3)
class(gsg)
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:",
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:",
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
dim(dataExpr)
head(dataExpr)
sampleTree = hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
powers = c(c(1:10), seq(from = 12, to=20, by=2))
type = 'unsigned' 
sft = pickSoftThreshold(dataExpr,powerVector= powers,
                        networkType=type, verbose=5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.9,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
     cex=cex1, col="red")
power = sft$powerEstimate
if (is.na(power)){
  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                 ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                        ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                               ifelse(type == "unsigned", 6, 12))       
                 )
  )
}
corType <- 'pearson'
type <- 'unsigned'
mingene <-10
cor <- WGCNA::cor
net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                       TOMType = type, minModuleSize = mingene,
                       reassignThreshold = 0, mergeCutHeight = 0.2,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = corType, 
                       loadTOMs=TRUE,
                       saveTOMFileBase = paste0('exprMat', ".tom"),
                       verbose = 3)
table(net$colors)
sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
table(mergedColors)
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.2,
                    addGuide = TRUE, guideHang = 0.05)
MEs = net$MEs
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)
trait <- 'clinical.csv'
if(trait != "") {
  traitData <- read.csv(file=trait, header=T, row.names=1,
  )
  sampleName = rownames(traitData)
  traitData = traitData[match(sampleName, rownames(traitData)), ]
}
MEs_colpheno = orderMEs(cbind(MEs_col, traitData))
plotEigengeneNetworks(MEs_colpheno, "Eigengene adjacency heatmap",
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T,
                      xLabelsAngle = 90)
TOM = TOMsimilarityFromExpr(dataExpr, power=power, corType=corType, networkType=type)
load(net$TOMFiles[1], verbose=T)
TOM <- as.matrix(TOM)
dissTOM = 1-TOM
plotTOM = dissTOM^7
diag(plotTOM) = NA
probes = colnames(dataExpr)
dimnames(TOM) <- list(probes, probes)
dim(dataExpr)
if (corType=="pearson") {
  modTraitCor = cor(MEs_col, traitData, use = "p")
  modTraitP = corPvalueStudent(modTraitCor, nSamples)
} else {
  modTraitCorP = bicorAndPvalue(MEs_col, traitData, robustY=robustY)
  modTraitCor = modTraitCorP$pearson
  modTraitP   = modTraitCorP$p
}
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 4), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
ZhongbingYang_color = colorRampPalette(c('#0899ba', '#FFFEFE','#e01e37'))(50)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitData), 
               yLabels = colnames(MEs_col), 
               cex.lab = 0.5, 
               ySymbols = colnames(MEs_col), colorLabels = FALSE, 
               colors = ZhongbingYang_color,
               textMatrix = textMatrix, setStdMargins = FALSE, 
               cex.text = 0.5, zlim = c(-0.2,0.2),
               main = paste("Module-trait relationships"))
table(net$colors)
sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
table(mergedColors)
library(dplyr)
module_trait_table <- data.frame(
  Module = rownames(modTraitCor),
  stringsAsFactors = FALSE
)
for (trait in colnames(modTraitCor)) {
  module_trait_table[[paste0(trait, "_r")]] <- modTraitCor[, trait]
  module_trait_table[[paste0(trait, "_p")]] <- modTraitP[, trait]
}
write.csv(
  module_trait_table,
  file = "ModuleTrait_Correlation_and_Pvalue.csv",
  row.names = FALSE
)
modNames = substring(names(MEs_col), 3)
if (corType=="pearson") {
  geneModuleMembership = as.data.frame(cor(dataExpr, MEs_col, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(
    as.matrix(geneModuleMembership), nSamples))
} else {
  geneModuleMembershipA = bicorAndPvalue(dataExpr, MEs_col, robustY=robustY)
  geneModuleMembership = geneModuleMembershipA$pearson
  MMPvalue   = geneModuleMembershipA$p
}
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
write.csv(MMPvalue,'3_MMPvalue.csv')
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]];
color<-unique(moduleColors)
for (i  in 1:length(color)) {
  y=t(assign(paste(color[i],"expr",sep = "."),dataExpr[moduleColors==color[i]]))
  write.csv(y,paste('4',color[i],"csv",sep = "."),quote = F)
}
traitNames=names(traitData)
geneTraitSignificance = as.data.frame(cor(dataExpr, traitData, use = "p", method="pearson")) 
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep="")
write.csv(GSPvalue,'5_GSPvalue.csv')
for (trait in traitNames){
  traitColumn=match(trait,traitNames)
  for (module in modNames){
    column = match(module, modNames)
    moduleGenes = moduleColors==module
    if (nrow(geneModuleMembership[moduleGenes,]) > 1){
      pdf(file=paste("9_", trait, "_", module,"_Module membership vs gene significance.pdf",sep=""),width=7,height=7)
      par(mfrow = c(1,1))
      verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                         abs(geneTraitSignificance[moduleGenes, traitColumn]),
                         xlab = paste("Module Membership in", module, "module"),
                         ylab = paste("Gene significance for ",trait),
                         main = paste("Module membership vs. gene significance\n"),
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
      dev.off()
    }
  }
}
names(dataExpr)
probes = names(dataExpr)
geneInfo0 = data.frame(probes= probes,
                       moduleColor = moduleColors)
for (Tra in 1:ncol(geneTraitSignificance))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneTraitSignificance[,Tra],
                         GSPvalue[, Tra])
  names(geneInfo0) = c(oldNames,names(geneTraitSignificance)[Tra],
                       names(GSPvalue)[Tra])
}
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[,mod],
                         MMPvalue[, mod])
  names(geneInfo0) = c(oldNames,names(geneModuleMembership)[mod],
                       names(MMPvalue)[mod])
}
geneOrder =order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]
write.table(geneInfo, file = "6_GS_and_MM.xls",sep="\t",row.names=F)
connectivity=abs(cor(dataExpr,use="p"))^power 
Alldegrees=intramodularConnectivity(connectivity, mergedColors) 
#Alldegrees$gene=rownames(Alldegrees)
dataKME=signedKME(dataExpr, MEs_col, outputColumnName="kME_MM.")
GS=cor(dataExpr, traitData, use="p") 
combin_data=cbind(Alldegrees,dataKME,GS)
write.csv(dataKME,"7_combine_hub-2.csv")
traitData <- traitData[rownames(dataExpr), , drop = FALSE]
nSamples <- nrow(dataExpr)
corType <- "pearson"
if (corType == "pearson") {
  modTraitCor <- cor(MEs_col, traitData, use = "pairwise.complete.obs")
  modTraitP   <- corPvalueStudent(modTraitCor, nSamples)
}
out_dir <- "WGCNA_Module_Trait_Scatter_PDF_CI"
dir.create(out_dir, showWarnings = FALSE)
module_names <- colnames(MEs_col)
trait_names  <- colnames(traitData)
for (module in module_names) {
  for (trait in trait_names) {
    x <- traitData[, trait]
    y <- MEs_col[, module]
    idx <- complete.cases(x, y)
    if (sum(idx) < 10) next
    x_use <- x[idx]
    y_use <- y[idx]
    r <- modTraitCor[module, trait]
    p <- modTraitP[module, trait]
    p_sci <- sprintf("%.4e", p)
    fit <- lm(y_use ~ x_use)
    pred <- predict(fit, interval = "confidence", level = 0.95)
    pdf_file <- file.path(
      out_dir,
      sprintf("%s_vs_%s_CI.pdf", module, trait)
    )
    pdf(pdf_file, width = 6, height = 6)
    plot(
      x_use, y_use,
      xlab = trait,
      ylab = module,
      main = paste(module, "~", trait),
      pch  = 19,
      col  = "steelblue",
      cex  = 1.2
    )
    abline(fit, col = "red", lwd = 2)
    polygon(
      c(x_use, rev(x_use)),
      c(pred[, "lwr"], rev(pred[, "upr"])),
      col = adjustcolor("gray", alpha.f = 0.35),
      border = NA
    )
    legend(
      "topright",
      legend = sprintf(
        "r = %.3f\np = %s\nn = %d",
        r, p_sci, sum(idx)
      ),
      bty = "n",
      cex = 1.1
    )
    dev.off()
  }
}
cor_table <- data.frame(
  Module = rep(rownames(modTraitCor), times = ncol(modTraitCor)),
  Trait  = rep(colnames(modTraitCor), each = nrow(modTraitCor)),
  r      = as.vector(modTraitCor),
  p      = as.vector(modTraitP)
)
cor_table$p_scientific <- sprintf("%.4e", cor_table$p)
cor_table$significance <- cut(
  cor_table$p,
  breaks = c(0, 0.001, 0.01, 0.05, 1),
  labels = c("***", "**", "*", "ns"),
  include.lowest = TRUE
)
write.csv(
  cor_table,
  "Module_Trait_Correlations_Pearson_WGCNA_aligned.csv",
  row.names = FALSE
)
salmon_module <- colnames(MEs_col)[grep("magenta", colnames(MEs_col))]
if (length(salmon_module) == 1) {
  mmse_col <- colnames(traitData)[grep("MMSE", colnames(traitData))]
  
  r     <- modTraitCor[salmon_module, mmse_col]
  p     <- modTraitP[salmon_module, mmse_col]
  p_sci <- sprintf("%.4e", p)
  
  cat("\n========== salmon 模块 vs MMSE ==========\n")
  cat("r     = ", r, "\n")
  cat("p     = ", p, "\n")
  cat("p_sci = ", p_sci, "\n")
  cat("=========================================\n\n")
}   
target_module <- "magenta"                    
stage_file    <- "disease_stage.csv"          
stage_order   <- c("Mi", "Mo", "S")         
ME_colname <- paste0("ME", target_module)
if (!ME_colname %in% colnames(MEs_col)) {
  stop(paste("找不到模块", ME_colname, "，可用的模块有：", paste(colnames(MEs_col), collapse=", ")))
}
ME_module <- MEs_col[[ME_colname]]
names(ME_module) <- rownames(MEs_col)
if (!file.exists(stage_file)) {
  stop(paste("找不到文件:", stage_file))
}
stage_info <- read.csv(stage_file, stringsAsFactors = FALSE)
if (!all(c("SampleID", "Stage") %in% colnames(stage_info))) {
  stop("分组文件必须包含 'SampleID' 和 'Stage' 两列")
}
common <- intersect(names(ME_module), stage_info$SampleID)
if (length(common) == 0) {
  stop("样本名不匹配，请检查")
}
ME_common    <- ME_module[common]
stage_common <- stage_info$Stage[match(common, stage_info$SampleID)]
plot_df <- data.frame(
  Sample = common,
  ME     = ME_common,
  Stage  = factor(stage_common, levels = stage_order)
)
plot_df <- na.omit(plot_df)
print(table(plot_df$Stage))
kwt <- kruskal.test(ME ~ Stage, data = plot_df)
cat("\n========== Kruskal-Wallis  ==========\n")
cat(sprintf("p-value = %.4e\n", kwt$p.value))
cat("=========================================\n\n")
if (!require(FSA)) install.packages("FSA")
library(FSA)
dunn_res <- dunnTest(ME ~ Stage, data = plot_df, method = "bh")
cat("========== Dunn （FDR）==========\n")
print(dunn_res$res)
cat("==================================================\n")
library(ggplot2)
library(ggpubr)
fill_colors <- setNames(c("magenta", "magenta", "magenta"), stage_order)
p <- ggplot(plot_df, aes(x = Stage, y = ME, fill = Stage)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6, shape = 21, color = "black", fill = "black") +
  stat_summary(fun = median, geom = "point", shape = 4, size = 3, color = "red") +
  labs(title = paste(target_module, "module eigengene across disease stages"),
       y = "Module eigengene (WGCNA)", x = "Disease stage") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_fill_manual(values = fill_colors)
sig_pairs <- dunn_res$res$Comparison[dunn_res$res$P.adj < 0.05]
if (length(sig_pairs) > 0) {
  comp_list <- lapply(strsplit(sig_pairs, " - "), function(x) list(x[1], x[2]))
  p <- p + stat_compare_means(comparisons = comp_list,
                              method = "wilcox.test",
                              label = "p.signif", hide.ns = TRUE,
                              tip.length = 0.01)
}
ggsave(paste0(target_module, "_ME_Boxplot.pdf"), plot = p, width = 5.5, height = 4.5)
ggsave(paste0(target_module, "_ME_Boxplot.png"), plot = p, width = 5.5, height = 4.5, dpi = 300)
cat("\n✅ 完成！图已保存：", paste0(target_module, "_ME_Boxplot.png"), "\n")
