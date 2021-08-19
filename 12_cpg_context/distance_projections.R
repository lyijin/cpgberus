
library(stringdist)
library(stringr)
library(RColorBrewer)
library(edgeR)
library(limma)


#####  Paths  #####

path_base = getwd()
path_cpgberus = "/media/hb-stopwatch/work/cpgberus"
path_txdb = file.path(path_base, "01_txdb")
path_cpgberus_covs = file.path(path_cpgberus, "04_parse_bismark_covs")
path_cpgberus_epic = file.path(path_cpgberus, "04_parse_methylationEPIC")


# TxDB
source(file.path(path_txdb, "annotate_ranges.R"))


"grch38p13.nnncgnnn_universe.tsv"



# library coverages
lib_covs = sapply(covs_grl, function(x) sum(x$meth_cov) + sum(x$umeth_cov))
lib_covs / 1e6

# Find all motifs across the samples
all_seqs = lapply(covs_grl, function(x) x$cpg_context_nnncgnnn)
#all_seqs = lapply(all_seqs, substr, 1, 3)
#all_seqs = lapply(all_seqs, substr, 5, 8)
seqs = unique(unlist(lapply(all_seqs, unique)))

# Fill a matrix of counts
seqs_mat = matrix(NA_real_,
                  nrow = length(seqs),
                  ncol = length(covs_grl),
                  dimnames = list(seqs, names(covs_grl)))

for(s in names(all_seqs)) {
  seqs_table = as.data.frame(table(all_seqs[[s]]))
  seqs_mat[as.character(seqs_table$Var1), s] = seqs_table$Freq
}

# Discard motifs with N in them
seqs_mat = seqs_mat[!grepl("N", rownames(seqs_mat)), ]

colSums(seqs_mat, na.rm = TRUE)

# Filter low counts
seqs_mat_sml = na.omit(seqs_mat)
#is_low = rowSums(seqs_mat, na.rm = TRUE) <= 20
#table(is_low)
#seqs_mat_sml = seqs_mat[!is_low, ]

log_seqs = log2(seqs_mat_sml)
plot(density(log_seqs[, 5]))

# Diff Exp, Limma

pheno = data.frame(sample = sub("(WR0\\d{2})V(\\d)(\\w)", "\\1", colnames(seqs_mat_sml)),
                   visit = as.integer(sub("WR0(\\d{2})V(\\d)(\\w)", "\\2", colnames(seqs_mat_sml))),
                   library = sub("WR0(\\d{2})V(\\d)(\\w)", "\\3", colnames(seqs_mat_sml)),
                   row.names = colnames(seqs_mat_sml))

mm = model.matrix(~ visit + library, data=pheno)

corfit <- duplicateCorrelation(log_seqs, mm, block=pheno$sample)
corfit$consensus
# Current value for model ~ 0 + Treat + gender + age_exact_yr + log(hsCRP_mg_l) + CD8T is: 

# design <- model.matrix(~ 0 + Treat + gender + Neutro + Mono, data=pd) # cor=0.3738861
# design <- model.matrix(~ 0 + Treat + gender + Neutro, data=pd) # cor=0.380095
fit1 <- lmFit(seqs_mat_sml, mm, block=pheno$sample, correlation=corfit$consensus)
fit1 <- eBayes(fit1)

#fit1$genes <- as.data.frame(annoEPIC[match(rownames(fit1), annoEPIC$Name), c(1:4,19,22:23, 35)])
dt1 <- decideTests(fit1, adjust.method="BH", p.value=0.01)
summary(dt1)
topTable(fit1, coef="libraryW", p.value=0.01, number = nrow(fit1), lfc = 0)

#seqs_mat_sml["GGCCGGGC", ]


covs_ol = covs_grl$WR025V1E[covs_grl$WR025V1E$cpg_context_nnncgnnn %in% c("GTCCGGGA", "TCCCGGAC")]

x = annotateMeth(covs_ol, tx = gencode_tx, cpgIslands = cpgIslands, cpgShores = cpgShores)


# MDS


seqs_dist = stringdistmatrix(rownames(log_seqs), method="hamming", useNames = TRUE)
seqs_mds = cmdscale(seqs_dist)


gc_content = factor(str_count(rownames(log_seqs), "G|C"))
gc_cols = brewer.pal(length(levels(gc_content)), "Spectral")[gc_content]
plot(seqs_mds, type="n")
text(seqs_mds, labels = rownames(seqs_mds), cex=0.7, col=gc_cols)

odds = seq(1, 8, by=2)
#diffs = rowMeans(log_seqs[, odds]) - rowMeans(log_seqs[, (odds + 1)])
#diffs = rowMeans(log_seqs[, odds] - log_seqs[, (odds + 1)])
diffs = rowMeans(log_seqs[, odds] / log_seqs[, (odds + 1)])
diffs = rowMeans(seqs_mat_sml[, odds] / seqs_mat_sml[, (odds + 1)])
#diffs = colMeans(seqs_mat_sml[, odds]) - colMeans(seqs_mat_sml[, (odds + 1)])
summary(diffs)
top_motifs = sort(diffs, decreasing = TRUE)[1:40]
top_motifs
seqs_mat_sml[names(top_motifs), ]
diff_cols = brewer.pal(7, "Greys")[cut(diffs, 7)]
plot(seqs_mds, type="n")
text(seqs_mds, labels = rownames(seqs_mds), cex=0.4, col=diff_cols)
#points(seqs_mds, cex=0.1, col=diff_cols)



seqs_mat_sml["GTCCGGGA", ]



seqs_clust = hclust(seqs_dist)
seqs_dend = as.dendrogram(seqs_clust)


# Normalise for library size
seqs_mat_norm = seqs_mat
for(i in 1:ncol(seqs_mat_norm)) {
  seqs_mat_norm[, i] = seqs_mat_norm[, i] / sum(seqs_mat_norm[, i], na.rm = TRUE) * 1e6
}

seqs_mat_norm2 = seqs_mat / colSums(seqs_mat, na.rm = TRUE) * 1e6
all.equal(seqs_mat_norm, seqs_mat_norm2)


heatmap(seqs_mat_norm2, Rowv = seqs_dend, Colv = NULL)

xx = dist(t(na.omit(seqs_mat_norm)))
seqs_mds_1 = cmdscale(xx)
plot(seqs_mds_1, type="n")
text(seqs_mds_1, labels = rownames(seqs_mds_1), cex=0.4)



# Diff Exp, voom
# Following, https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
d0 <- DGEList(na.omit(seqs_mat))
d0 <- calcNormFactors(d0)
d0

cutoff <- 10
drop <- which(apply(cpm(d0), 1, max) < cutoff)
table(drop)
d <- d0[-drop,] 
dim(d) 
plotMDS(d)


y <- voom(d, mm, plot = T)

# Doesn't look like there's a mean-variance trend we need to model away with voom



