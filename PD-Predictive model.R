rm(list = ls())
if (dev.cur() > 1) dev.off()
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(reshape2)

input_file <- "1.csv"
k_fold   <- 5          
top_n    <- 30         
ntree_rf <- 1000  
seed_split <- 100
seed_full  <- 100
seed_cv    <- 100
seed_glm   <- 100
fig_dir <- "figures"
if (!dir.exists(fig_dir)) dir.create(fig_dir)
impute_zero <- function(df, features) {
  for (v in features) {
    x <- df[[v]]
    x[is.na(x)] <- 0
    df[[v]] <- x
  }
  df
}
plot_confusion_matrix <- function(actual, predicted, title = "",
                                  file_name = "confusion_matrix",
                                  save_dir  = fig_dir) {
  actual    <- factor(actual,    levels = c("Ctrl", "PD"))
  predicted <- factor(predicted, levels = c("Ctrl", "PD"))
  cm <- caret::confusionMatrix(predicted, actual, positive = "PD")
  sens <- round(as.numeric(cm$byClass["Sensitivity"]), 4)
  spec <- round(as.numeric(cm$byClass["Specificity"]), 4)
  ppv  <- round(as.numeric(cm$byClass["Pos Pred Value"]), 4)
  npv  <- round(as.numeric(cm$byClass["Neg Pred Value"]), 4)
  acc  <- round(as.numeric(cm$overall["Accuracy"]), 4)
  tab <- as.data.frame(cm$table)
  colnames(tab) <- c("Prediction", "Reference", "Freq")
  tab$Cell <- with(tab,
                   ifelse(Prediction == "PD"   & Reference == "PD",   "TP",
                          ifelse(Prediction == "Ctrl" & Reference == "PD",   "FN",
                                 ifelse(Prediction == "PD"   & Reference == "Ctrl", "FP", "TN"))))
  tab$Label <- sprintf("%s\n%d", tab$Cell, tab$Freq)
  tab$Prediction <- factor(tab$Prediction, levels = c("Ctrl", "PD"))
  tab$Reference  <- factor(tab$Reference,  levels = c("Ctrl", "PD"))
  p <- ggplot(tab, aes(x = Reference, y = Prediction, fill = Freq)) +
    geom_tile(color = "white", linewidth = 1.2) +
    geom_text(aes(label = Label), size = 6, fontface = "bold", color = "grey15") +
    scale_fill_gradient(low = "#E8F4FA", high = "#4DA8DA") +
    labs(title = title,
         subtitle = sprintf("Sensitivity = %.4f | Specificity = %.4f\nPPV = %.4f | NPV = %.4f | Accuracy = %.4f",
                            sens, spec, ppv, npv, acc),
         x = "Actual (Reference)", y = "Predicted") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "none", panel.grid = element_blank())
  ggsave(file.path(save_dir, paste0(file_name, ".png")), plot = p,
         width = 6.5, height = 5.5, dpi = 300, bg = "white")
  ggsave(file.path(save_dir, paste0(file_name, ".pdf")), plot = p,
         width = 6.5, height = 5.5, device = cairo_pdf, bg = "white")
  print(p)
  cm
}
plot_roc_curve <- function(gold, prob, title = "",
                           file_name = "roc_curve",
                           save_dir  = fig_dir) {
  gold    <- factor(gold, levels = c("Ctrl", "PD"))
  roc_obj <- pROC::roc(gold, prob, levels = c("Ctrl", "PD"),
                       direction = "<", quiet = TRUE, ci = TRUE)
  draw_fun <- function() {
    plot(1 - roc_obj$specificities, roc_obj$sensitivities,
         col = "#0093FF", lwd = 2.5, type = "l", xlim = c(0, 1), ylim = c(0, 1),
         xlab = "1 - Specificity (FPR)", ylab = "Sensitivity (TPR)", main = title)
    abline(0, 1, lty = 2, lwd = 1.5, col = "grey40")
    legend("bottomright", bty = "n",
           legend = sprintf("AUC: %.4f (%.4f-%.4f)",
                            as.numeric(roc_obj$auc), roc_obj$ci[1], roc_obj$ci[3]), cex = 1.0)
  }
  png(file.path(save_dir, paste0(file_name, ".png")),
      width = 6.5, height = 5.5, units = "in", res = 300); draw_fun(); dev.off()
  pdf(file.path(save_dir, paste0(file_name, ".pdf")), width = 6.5, height = 5.5)
  draw_fun(); dev.off()
  draw_fun()
  roc_obj
}
set.seed(123)
dat <- read.csv(input_file, header = TRUE, check.names = FALSE,
                row.names = 1, fileEncoding = "UTF-8-BOM")
if (!"Group" %in% names(dat)) stop("数据中找不到 Group 列")
dat$Group <- factor(dat$Group, levels = c("Ctrl", "PD"))
colnames(dat) <- gsub("-", "_", colnames(dat))
feat_all <- setdiff(names(dat), "Group")
dat[feat_all] <- lapply(dat[feat_all],
                        function(x) suppressWarnings(as.numeric(as.character(x))))
all_na <- sapply(dat[feat_all], function(x) all(is.na(x)))
dat <- dat[, c("Group", feat_all[!all_na]), drop = FALSE]
cat("样本数:", nrow(dat), "；特征数:", ncol(dat) - 1, "\n")
print(table(dat$Group))
set.seed(seed_split)
train_idx <- createDataPartition(dat$Group, p = 0.70, list = FALSE)
train <- dat[train_idx, , drop = FALSE]
test  <- dat[-train_idx, , drop = FALSE]
write.csv(train, "train70_raw.csv", row.names = TRUE)
write.csv(test,  "test30_raw.csv",  row.names = TRUE)
cat("训练集:", nrow(train), "；测试集:", nrow(test), "\n")
print(table(train$Group))
print(table(test$Group))
features <- setdiff(names(train), "Group")
train_zero <- impute_zero(train, features)
keep_full <- features[sapply(train_zero[features], function(x) length(unique(x)) > 1)]
set.seed(seed_full)
rf_full <- randomForest(x = train_zero[, keep_full, drop = FALSE],
                        y = train_zero$Group, ntree = ntree_rf,
                        importance = TRUE, proximity = FALSE)
imp_full <- importance(rf_full, type = 2)
imp_full <- imp_full[order(imp_full[, "MeanDecreaseGini"], decreasing = TRUE), , drop = FALSE]
cat("\n===== 【RF 特征排序】全训练集 MeanDecreaseGini 前 30 =====\n")
print(head(imp_full, 30))
write.csv(imp_full, "RF_importance_MeanDecreaseGini.csv", row.names = TRUE)
png(file.path(fig_dir, "RF_varImpPlot.png"), width = 900, height = 700)
varImpPlot(rf_full, main = "Random Forest Variable Importance (Feature Ranking Only)")
dev.off()
mtry_rf <- max(1, floor(sqrt(length(features))))
set.seed(seed_cv)
folds <- createFolds(train$Group, k = k_fold, list = TRUE, returnTrain = FALSE)
gini_list <- list()
for (k in seq_along(folds)) {
  cat("正在运行 CV fold:", k, "\n")
  cv_train <- impute_zero(train[-folds[[k]], , drop = FALSE], features)
  keep <- features[sapply(cv_train[features], function(x) length(unique(x)) > 1)]
  x_cv <- cv_train[, keep, drop = FALSE]
  y_cv <- cv_train$Group
  set.seed(1000 + k)
  rf_cv <- randomForest(x = x_cv, y = y_cv, ntree = ntree_rf,
                        mtry = min(mtry_rf, ncol(x_cv)), importance = TRUE)
  imp <- importance(rf_cv, type = 2)[, "MeanDecreaseGini"]
  gini_full <- setNames(rep(0, length(features)), features)
  gini_full[names(imp)] <- imp
  gini_list[[paste0("Fold", k)]] <- gini_full
}
gini_mat  <- do.call(cbind, gini_list)
fold_cols <- colnames(gini_mat)
gini_df <- data.frame(Feature = rownames(gini_mat), gini_mat, check.names = FALSE)
gini_df$Mean_Gini <- rowMeans(gini_df[, fold_cols, drop = FALSE])
gini_df$SD_Gini   <- apply(gini_df[, fold_cols, drop = FALSE], 1, sd)
rank_mat <- apply(gini_df[, fold_cols, drop = FALSE], 2,
                  function(x) rank(-x, ties.method = "min"))
gini_df$Mean_Rank <- rowMeans(rank_mat)
gini_df$TopN_Freq <- rowSums(rank_mat <= top_n)
gini_df <- gini_df[order(-gini_df$Mean_Gini), ]
write.csv(gini_df, "gini_cv_summary.csv", row.names = FALSE)
gini_long <- melt(gini_df[, c("Feature", fold_cols), drop = FALSE],
                  id.vars = "Feature", variable.name = "Fold",
                  value.name = "MeanDecreaseGini")
write.csv(gini_long, "gini_cv_long.csv", row.names = FALSE)
top_plot_n <- min(30, nrow(gini_df))
p_bar <- ggplot(gini_df[1:top_plot_n, ],
                aes(x = reorder(Feature, Mean_Gini), y = Mean_Gini)) +
  geom_col(fill = "#0093FF") + coord_flip() + theme_bw() +
  labs(x = "", y = "Mean Decrease Gini") +
  ggtitle("RF 5-fold CV MeanDecreaseGini (feature ranking, zero imputation)")
ggsave(file.path(fig_dir, "gini_mean_top30.pdf"), p_bar, width = 7, height = 6)
plot_df <- gini_long[gini_long$Feature %in% gini_df$Feature[1:top_plot_n], ]
p_heat <- ggplot(plot_df, aes(x = Feature, y = Fold, fill = MeanDecreaseGini)) +
  geom_tile() + scale_fill_gradient(low = "white", high = "red") +
  theme_bw() + theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(x = "", y = "", fill = "Gini")
ggsave(file.path(fig_dir, "gini_cv_heatmap_top30.pdf"), p_heat, width = 10, height = 4)
manual_features <- c("BLM")
feature_source  <- "manual_specification"
manual_features <- gsub("-", "_", manual_features)
if (length(manual_features) == 0) {
  stop(
    "\n================ 已输出 Gini 表，等待你填写特征 ================\n",
    "请打开以下文件查看 Gini 排名：\n",
    "  RF_importance_MeanDecreaseGini.csv\n",
    "  gini_cv_summary.csv\n",
    "  figures/gini_mean_top30.pdf\n",
    "  figures/gini_cv_heatmap_top30.pdf\n\n",
    "然后在代码中部的【人工决策区】填写：\n",
    "  manual_features <- c(\"蛋白1\", \"蛋白2\", \"蛋白3\")\n",
    "再重新运行本脚本，才会进入逻辑回归建模。\n",
    "==============================================================\n"
  )
}
missing_feat <- setdiff(manual_features, colnames(dat))
if (length(missing_feat) > 0) {
  stop("以下手动指定特征在数据中不存在，请检查拼写：\n",
       paste(missing_feat, collapse = ", "))
}
selected <- manual_features
cat("\n===== 最终进入逻辑回归的蛋白 =====\n")
print(selected)
cat("特征数:", length(selected), "\n")
cat("特征来源:", feature_source, "\n")
write.csv(data.frame(Feature = selected, Source = feature_source),
          "manual_selected_features.csv", row.names = FALSE)
train_final <- impute_zero(train, selected)
test_final  <- impute_zero(test,  selected)
keep_final <- selected[sapply(train_final[selected], function(x) length(unique(x)) > 1)]
if (length(keep_final) == 0) stop("手动指定特征在训练集中均为零方差")
if (length(keep_final) < length(selected)) {
  warning("以下特征在训练集中零方差，已被剔除：\n",
          paste(setdiff(selected, keep_final), collapse = ", "))
}
selected <- keep_final
cat("\n最终建模特征数:", length(selected), "\n")
print(selected)
write.csv(data.frame(Feature = selected), "manual_selected_features_final.csv", row.names = FALSE)
form_final <- as.formula(
  paste("Group ~", paste(sprintf("`%s`", selected), collapse = " + "))
)
cat("\n逻辑回归公式：", deparse(form_final), "\n")
cat("\n===== 入选特征在两组中的 0 比例 =====\n")
zero_ratio <- sapply(selected,
                     function(v) tapply(train_final[[v]] == 0, train_final$Group, mean))
print(round(zero_ratio, 3))
fitControl <- trainControl(method = "repeatedcv", number = 10, repeats = 10,
                           classProbs = TRUE, savePredictions = TRUE,
                           summaryFunction = twoClassSummary)
set.seed(seed_glm)
glm.model.cv <- train(form_final, data = train_final, method = "glm",
                      family = "binomial", metric = "ROC", trControl = fitControl)
print(glm.model.cv)
cat("\n最终逻辑回归模型摘要：\n")
print(summary(glm.model.cv$finalModel))
sink("final_glm_summary.txt"); print(summary(glm.model.cv$finalModel)); sink()
fitted_train <- fitted(glm.model.cv$finalModel)
if (any(fitted_train > 1 - 1e-8 | fitted_train < 1e-8)) {
  cat("\n[警告] 训练集拟合概率接近 0 或 1，可能存在完全分离。\n")
} else {
  cat("\n未检测到明显完全分离。\n")
}
glm.probs_train <- predict(glm.model.cv, train_final, type = "prob")[, "PD"]
glm.pred_train  <- predict(glm.model.cv, train_final)
cm_train <- plot_confusion_matrix(train_final$Group, glm.pred_train,
                                  "Confusion Matrix - Training Set", file_name = "01_CM_train")
print(cm_train)
roc_train <- plot_roc_curve(train_final$Group, glm.probs_train,
                            "ROC - Training Set", file_name = "02_ROC_train")
cat("训练集表观 AUC:", round(as.numeric(roc_train$auc), 4), "\n")
write.csv(data.frame(Sample = rownames(train_final), Prob = glm.probs_train,
                     Group = train_final$Group, Pred = glm.pred_train),
          "train_pred_prob.csv", row.names = FALSE)
glm.probs_test <- predict(glm.model.cv, test_final, type = "prob")[, "PD"]
glm.pred_test  <- predict(glm.model.cv, test_final)
cm_test <- plot_confusion_matrix(test_final$Group, glm.pred_test,
                                 "Confusion Matrix - Test Set", file_name = "03_CM_test")
print(cm_test)
roc_test <- plot_roc_curve(test_final$Group, glm.probs_test,
                           "ROC - Test Set", file_name = "04_ROC_test")
cat("测试集 AUC:", round(as.numeric(roc_test$auc), 4), "\n")
write.csv(data.frame(Sample = rownames(test_final), Prob = glm.probs_test,
                     Group = test_final$Group, Pred = glm.pred_test),
          "test_pred_prob.csv", row.names = FALSE)
write.csv(train_final, "train_split_zero.csv", row.names = TRUE)
write.csv(test_final,  "test_split_zero.csv",  row.names = TRUE)
cv_results <- data.frame(CV_AUC_Mean = round(glm.model.cv$results$ROC, 4),
                         CV_AUC_SD   = round(glm.model.cv$results$ROCSD, 4))
cat("\n===== 逻辑回归 10 折重复 CV 最终结果 =====\n"); print(cv_results)
write.csv(cv_results, "cv_results.csv", row.names = FALSE)
write.csv(glm.model.cv$resample, "cv_resample_detail.csv", row.names = FALSE)
write.csv(varImp(glm.model.cv)$importance, "glm_varImp.csv", row.names = TRUE)
auc_summary <- data.frame(
  Dataset = c("Train", "Test"), Model = "glm_logistic",
  AUC = c(as.numeric(roc_train$auc), as.numeric(roc_test$auc)),
  CI_low  = c(roc_train$ci[1], roc_test$ci[1]),
  CI_high = c(roc_train$ci[3], roc_test$ci[3]))
write.csv(auc_summary, "AUC_summary.csv", row.names = FALSE)
while (dev.cur() > 1) dev.off()
cat("\n================ 运行完成 ================\n")
cat("RF 仅用于特征排序；逻辑回归为最终分类器。\n")
cat("本次建模实际使用特征: manual_selected_features_final.csv\n")
cat("Gini 参考表: RF_importance_MeanDecreaseGini.csv / gini_cv_summary.csv\n")
cat("CV 结果: cv_results.csv / cv_resample_detail.csv\n")
cat("AUC 汇总: AUC_summary.csv\n")
cat("图形目录:", normalizePath(fig_dir), "\n")

