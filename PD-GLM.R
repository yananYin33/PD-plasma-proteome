
rm(list = ls())
try(dev.off(), silent = TRUE)  
library(ggpubr)
data1 = read.csv("1.csv", header = T, row.names = 1)
OUTDIR <- "."                  ## ← 输出目录，改成你需要的路径（如 "/data/workspace"）
YVAR <- "HY"         
XVAR <- "APLP1"   
dat <- data1
names(dat) <- c("year", "HY", "APLP1")
dat$yr        <- as.integer(sub("^Y", "", dat$year))          # Y1 -> 1
dat$year      <- factor(dat$year,
                        levels = paste0("Y", sort(unique(dat$yr))))
dat$HY        <- as.numeric(dat$HY)
dat$APLP1 <- as.numeric(dat$APLP1)
n0 <- nrow(dat)
dat <- dat[complete.cases(dat[, c("year", "yr", "HY", "APLP1")]), ]
cat("总样本量 N =", nrow(dat), " （剔除不完整记录", n0 - nrow(dat), "条）\n\n")
cat("各年份样本量：\n"); print(table(dat$year)); cat("\n")
cat("HY 分期分布：\n"); print(table(dat$HY)); cat("\n")
F_BASE <- as.formula(paste(YVAR, "~", XVAR))
F_MAIN <- as.formula(paste(YVAR, "~", XVAR, "+ yr"))   
F_INT  <- as.formula(paste(YVAR, "~", XVAR, "* yr"))   
fit_one <- function(s, formula = F_BASE) {
  g  <- glm(formula, data = s, family = gaussian(link = "identity"))
  b  <- coef(summary(g))
  nm <- XVAR
  ci <- confint.default(g)[nm, ]
  # ====
  r  <- cor(s[[XVAR]], s[[YVAR]], use = "complete.obs")
  r2 <- 1 - g$deviance / g$null.deviance         
  data.frame(beta = b[nm, "Estimate"],
             se   = b[nm, "Std. Error"],
             lo95 = ci[1],
             hi95 = ci[2],
             t    = b[nm, "t value"],
             p    = b[nm, "Pr(>|t|)"],
             r    = r,                            
             r2   = r2)                         
}
res <- do.call(rbind, lapply(levels(dat$year), function(y) {
  s <- dat[dat$year == y, ]
  cbind(data.frame(year = y, n = nrow(s)), fit_one(s))
}))
res$p_bonf <- p.adjust(res$p, method = "bonferroni")  
res$p_fdr  <- p.adjust(res$p, method = "fdr")
cat("================================================================\n")
cat("主分析：逐年 GLM（gaussian / identity）  ", deparse(F_BASE), "\n")
cat("================================================================\n")
print(format(res, digits = 4))
cat("\n")
write.csv(res, file.path(OUTDIR, "GLM_per_year.csv"), row.names = FALSE)
cat("================================================================\n")
cat("等价性自检：glm(gaussian/identity)  vs  lm()\n")
cat("================================================================\n")
chk <- do.call(rbind, lapply(levels(dat$year), function(y) {
  s <- dat[dat$year == y, ]
  g <- glm(F_BASE, data = s, family = gaussian(link = "identity"))
  o <- lm(F_BASE,  data = s)
  data.frame(year   = y,
             d_beta = abs(coef(g)[XVAR] - coef(o)[XVAR]),
             d_se   = abs(coef(summary(g))[XVAR, 2] - coef(summary(o))[XVAR, 2]),
             d_p    = abs(coef(summary(g))[XVAR, 4] - coef(summary(o))[XVAR, 4]),
             d_ci   = max(abs(confint.default(g)[XVAR, ] - confint(o)[XVAR, ])))
}))
print(format(chk, digits = 3))
cat("\n>>> 各列应全为 0（浮点误差 < 1e-10），即 GLM 与原有 OLS 结果完全一致\n")
cat(">>> 若不为 0，请检查是否误用了 confint()（profile）而非 confint.default()\n\n")
cat("================================================================\n")
cat("敏感性分析：合并模型（全部年份）\n")
cat("================================================================\n")
m_main <- glm(F_MAIN, data = dat, family = gaussian(link = "identity"))
m_int  <- glm(F_INT,  data = dat, family = gaussian(link = "identity"))
# =====
r2_main <- 1 - m_main$deviance / m_main$null.deviance
r2_int  <- 1 - m_int$deviance  / m_int$null.deviance
cat("--- 模型 1：", deparse(F_MAIN), "---\n")
print(round(coef(summary(m_main)), 4))
cat("R² =", round(r2_main, 4), "\n\n")
cat("--- 模型 2：", deparse(F_INT), "---\n")
print(round(coef(summary(m_int)), 4))
cat("R² =", round(r2_int, 4), "\n\n")
cat("--- 似然比检验：交互项是否必要 ---\n")
print(anova(m_main, m_int, test = "Chisq"))
cat("\nAIC：主效应模型 =", round(AIC(m_main), 2),
    "  交互模型 =", round(AIC(m_int), 2), "\n\n")
write.csv(round(coef(summary(m_int)), 4),
          file.path(OUTDIR, "GLM_pooled_interaction.csv"))
cat("================================================================\n")
cat("可选敏感性分析：有序 logistic（累积 logit，OR 口径）\n")
cat("================================================================\n")
if (requireNamespace("MASS", quietly = TRUE)) {
  ord <- do.call(rbind, lapply(levels(dat$year), function(y) {
    s     <- dat[dat$year == y, ]
    s$hy  <- factor(s[[YVAR]], ordered = TRUE)
    if (nlevels(s$hy) < 2) return(NULL)
    ok <- try(suppressWarnings(
      MASS::polr(as.formula(paste("hy ~", XVAR)), data = s, Hess = TRUE)),
      silent = TRUE)
    if (inherits(ok, "try-error")) return(NULL)
    ct <- coef(summary(ok))
    data.frame(year = y, n = nrow(s),
               logOR = ct[XVAR, "Value"],
               OR    = exp(ct[XVAR, "Value"]),
               p     = ct[XVAR, "Pr(>|t|)"])
  }))
  print(format(ord, digits = 4))
  write.csv(ord, file.path(OUTDIR, "GLM_ordinal_sensitivity.csv"), row.names = FALSE)
} else {
  cat("MASS 不可用，已跳过\n")
}
cat("\n")
cat("================================================================\n")
cat("诊断：各年最大 Cook 距离（经验阈值 4/n）\n")
cat("================================================================\n")
ck <- do.call(rbind, lapply(levels(dat$year), function(y) {
  s  <- dat[dat$year == y, ]
  cd <- cooks.distance(lm(F_BASE, data = s))
  data.frame(year = y, n = nrow(s),
             max_cook  = round(max(cd), 3),
             threshold = round(4 / nrow(s), 3),
             flag      = ifelse(max(cd) > 4 / nrow(s), "需关注", ""),
             row.names = NULL)
}))
print(ck, row.names = FALSE)
cat("\n")
#=====
res$lab <- with(res, sprintf(
  "β = %.3f [%.3f, %.3f]\nP = %s\nR² = %.3f",
  beta, lo95, hi95,
  ifelse(p < .001, "< .001", sprintf("%.3f", p)),
  r2
))
labdat <- data.frame(
  year  = res$year,
  xpos  = max(dat[[XVAR]], na.rm = TRUE),
  ypos  = max(dat[[YVAR]], na.rm = TRUE),
  lab   = res$lab,
  n     = res$n
)
labdat$lab2 <- paste0("n = ", labdat$n, "\n", labdat$lab)
p <- ggscatter(
  dat, x = XVAR, y = YVAR,
  facet.by    = "year",
  ncol        = 4, nrow = 2,
  add         = "reg.line",
  conf.int    = TRUE,
  add.params  = list(color = "red", fill = "#d6e2f5", size = 1.5),
  color       = "black",
  alpha       = 0.85,
  size        = 1,
  xlab        = paste(XVAR, "at baseline"),
  ylab        = YVAR,
  ggtheme     = theme_bw()
) +
  geom_text(data = labdat,
            aes(x = xpos, y = ypos, label = lab2),
            inherit.aes = FALSE,
            hjust = 1.05, vjust = 1.3, size = 2.8, lineheight = 0.9) +
  theme(strip.background = element_rect(fill = "#f2f2f2", colour = "grey70"),
        strip.text       = element_text(face = "bold", size = 10),
        panel.grid.minor = element_blank())
print(p)
ggsave(file.path(OUTDIR, "Fig_GLM_panels.pdf"), p, width = 12, height = 6)
ggsave(file.path(OUTDIR, "Fig_GLM_panels.png"), p, width = 12, height = 6, dpi = 300)
cat("----------------------------------------------------------------\n")
cat("已输出：\n")
cat("  结果表   ", file.path(OUTDIR, "GLM_per_year.csv"), "\n")
cat("  合并模型 ", file.path(OUTDIR, "GLM_pooled_interaction.csv"), "\n")
cat("  有序模型 ", file.path(OUTDIR, "GLM_ordinal_sensitivity.csv"), "\n")
cat("  图形     ", file.path(OUTDIR, "Fig_GLM_panels.pdf / .png"), "\n")