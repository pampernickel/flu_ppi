##========================================================
##
##  Credits:
##  Code written by MPDobay for flu interaction network analysis
##  SystemsX project: BCF/SIB
##
##========================================================

# ---
# ---
# comparison of v9_05 and v10 networks based on host factors:
# Figure 2D
library(ggplot2)
library(igraph)
library(reshape2)

# Files created with ./paper.scripts/string_annotation.r
load("./r.data.files/string.v10.rda")
load("./r.data.files/string.v9.rda")
g -> g.10
g -> g.9

# get graph intersection
graph.intersection(g.10, g.9, keep.all.vertices = F) -> intersection
c("neighborhood", "cooccurence", "coexpression", "experiments", "experimental", "database",
  "textmining", "combined_score") -> names
get.data.frame(intersection, what="edges") -> g.df

unlist(sapply(names, function(x) grep(x, colnames(g.df)))) -> col.ind
g.df[,col.ind] -> g.df
colnames(g.df)[which(colnames(g.df) %in% "experiments")] <- "experiments_1"
colnames(g.df)[which(colnames(g.df) %in% "experimental")] <- "experiments_2"
c(grep("_1", colnames(g.df)),grep("_2", colnames(g.df))) -> col.ind
g.df[,col.ind] -> g.df
gsub("_1", ".v_10", colnames(g.df)) -> colnames(g.df)
gsub("_2", ".v_9_05", colnames(g.df)) -> colnames(g.df)

c("neighborhood", "cooccurence", "coexpression", "experiments", "database",
  "textmining", "combined_score") -> names

# put columns with the same name side by side
df.cor <- matrix(NA, nrow=0, ncol=3)
colnames(df.cor) <- c("v_10", "v_9_05", "variable")
for (i in 1:length(names)){
  cbind(g.df[,grep(names[i], colnames(g.df))],
        rep(names[i], nrow(g.df))) -> t
  colnames(t) <- c("v_10", "v_9_05", "variable")
  rbind(df.cor, t) -> df.cor
}

as.data.frame(df.cor) -> df.cor
as.numeric(as.character(df.cor[,1])) -> df.cor[,1]
as.numeric(as.character(df.cor[,2])) -> df.cor[,2]

colnames(df.cor) <- c("v10", "v9", "variable")

# subset, as this would otherwise the visualization would take too long
# and would result in a very large figure
sample(nrow(df.cor), 1437144, replace=F) -> ind #33% of total data points
ggplot(df.cor[ind,], aes(x=v10, y=v9))+
  geom_point(size=0.2)+facet_wrap(~variable)+
  theme(axis.text=element_text(size=12),
        strip.text.x=element_text(size=16),
        axis.title=element_text(size=16))+
  xlab("Version 10")+
  ylab("Version 9_05")