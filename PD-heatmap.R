library(readxl)
library(dplyr)
library(tidyr)
data <- read.csv('1.csv',row.names = 1,header = T)
bk = unique(c(seq(-1,1, length=50)))
library(pheatmap)
pdf("1.pdf",width = 5,height = 7)
p1 = pheatmap(data, 
              cluster_cols = F,
              breaks = bk,
              cluster_row = F,
              show_colnames  = T,
              show_rownames  = T,
              color = colorRampPalette(c("#2A539A",
                                         'white',
                                         '#D9362C'
                                         ))(50), 
              legend = T,
              annotation_legend = T,
              border_color = 'NA',
              treeheight_row=5,
              treeheight_col=5,
              fontsize = 10,
              fontsize_row=6
            )
dev.off()