
library(stringdist)
library(stringr)
library(GenomicRanges)
library(Biostrings)
library(RColorBrewer)
library(limma)


#####  Paths  #####

path_base = getwd()
# Assumes your ./data symlink points to hb-stopwatch/work/cpgberus
path_cpgberus = file.path(path_base, "data")
path_txdb = file.path(path_cpgberus, "01_txdb")
path_cpgberus_covs = file.path(path_cpgberus, "04_parse_bismark_covs")


#####  Load  #####

load(file=file.path(path_cpgberus_covs, "motif_matrices.RData"))


#####  Log transform  #####

# Log2 transformation reduces the skewness of the data somewhat
plot(density(motifs$all_enrich[, 1]))
plot(density(log2(motifs$all_enrich[, 1])))

# So let's transform it
enrich_names = names(motifs)[grep("_enrich", names(motifs))]
log_motifs = lapply(enrich_names, function(n) log2(motifs[[n]]))
names(log_motifs) = enrich_names
rm(enrich_names)


#####  Diff Exp  #####

# Make model
pheno = data.frame(sample = sub("(WR0\\d{2})V(\\d)(\\w)", "\\1", colnames(log_motifs$all_enrich)),
                   visit = as.integer(sub("WR0(\\d{2})V(\\d)(\\w)", "\\2", colnames(log_motifs$all_enrich))),
                   library = sub("WR0(\\d{2})V(\\d)(\\w)", "\\3", colnames(log_motifs$all_enrich)),
                   row.names = colnames(log_motifs$all_enrich))

# Differential expression conditional upon the visit number and library type
mm = model.matrix(~ visit + library, data=pheno)

# However, the data has some correlation, as the same person has multiple datapoints
# We need to address this correlation with a blocking variable in the design
corfit <- duplicateCorrelation(log_motifs$all_enrich, mm, block=pheno$sample)
corfit$consensus  # -0.1105491
fit1 <- lmFit(log_motifs$all_enrich, mm, block=pheno$sample, correlation=corfit$consensus)
fit1 <- eBayes(fit1)

# Which motifs are different?
dt1 <- decideTests(fit1, adjust.method="BH", p.value=0.01)
summary(dt1)
# Now a topTable for library type
# W (WGBS) is the standard case, so logFC < 0 means that motif was less in WGBS
# A logFC > 0 means the motif was greater than in WGBS than in EM-seq
tt = topTable(fit1, coef="libraryW", p.value=1e-5, number = nrow(fit1))
dim(tt)
tt

# Dom, it would be great to map your motif plots onto motifs that are higher/lower in some form of Limma toptable
# We might need to do this in combination with some form of clustering analysis. Either HClust and cutree or MDS for ordering?
# Alternatively, we can just tighten the screws on the p-value cutoff on the toptable

# Is up in EM-Seq. This confirms our hypothesis, that GC-rich motifs are under-represented in WGBS.
tt[tt$logFC < 0, ]

# Is up in WGBS. What does this mean? The TET enzyme doesn't like to bind here? Dunno...
tt[tt$logFC > 0, ]

# One big gotta we need to address which Yi Jin highlighted is the appearance of reverse complements in the table
dna = DNAStringSet(rownames(tt))
rc = match(dna, reverseComplement(dna))
# e.g. row 11 is the reverseComplement of row 1, row 9 of row 2, etc.
tt[c(1, 11), ]
tt[c(2, 9), ]

# We need to remove/combine all revcomps before making motif plots


#####  Log Ratios  #####

# One really useful property of logs is for difference ratios
# Say for one feature we have half in group A vs B and double for another feature
# By expressing as a log, things that are smaller or larger are symmetrical around 0
log2(1/2)
log2(2)

# The gotcha is that log doesn't like 0.
# The common stats hack for data that contains 0 is to add 0.5 or something like that to all the data


odds = seq(1, 8, by=2)
# EM-seq / WGBS
diffs = rowMeans(log_motifs$all_enrich[, odds] / log_motifs$all_enrich[, (odds + 1)])
rm(odds)
summary(diffs)
diffs = diffs[order(diffs, decreasing = TRUE)]

# There's more motifs enriched with EM-seq than WGBS
plot(diffs, type="l")
abline(h=c(0.985, 1, 1.015), col=c("orange", "red", "orange"))


wgbs_high = DNAStringSet(names(diffs)[diffs < 0.98])
emseq_high = DNAStringSet(names(diffs)[diffs > 1.02])


#####  MDS of the coverage data  #####

# Let's explore the data via MDS and look for global patterns

# By sample
seqs_mds = cmdscale(dist(t(log_motifs$all_enrich)))
# They split across dim1 by library type
# Pleasingly, the samples also split across dim2 by person id and visit
plot(seqs_mds, type="n")
text(seqs_mds, labels = rownames(seqs_mds), cex=0.8)

# By motif
seqs_mds = cmdscale(dist(log_motifs$all_enrich))
plot(seqs_mds, type="n")
text(seqs_mds, labels = rownames(seqs_mds), cex=0.4, col="#80808080")

# Make some vectors for colouring by beta or GC
mean_betas = rowMeans(motifs$beta[rownames(seqs_mds), ])
mean_gcs = rowMeans(motifs$gc_pct[rownames(seqs_mds), ])
beta_cols = brewer.pal(7, "Spectral")[cut(mean_betas, 7)]
gc_cols = brewer.pal(7, "Spectral")[cut(mean_gcs, 7)]

# Dim 1 & 2 is influenced by beta/GC content of the motifs
plot(seqs_mds, cex=0.5, col=beta_cols, main="Meth (%)", pch=16)
legend("topright", legend=levels(cut(mean_betas, 7)),
       fill=brewer.pal(7, "Spectral"))

plot(seqs_mds, cex=0.5, col=gc_cols, main="GC")
legend("topright", legend=levels(cut(mean_gcs, 7)),
       fill=brewer.pal(7, "Spectral"))

# Dim 1 & 2 is influenced by motif coverage
mean_log_cov = rowMeans(log_motifs$all_enrich)
cov_cols = brewer.pal(7, "Blues")[cut(mean_log_cov, 7)]
plot(seqs_mds, cex=0.5, col=cov_cols, main="Log2 Coverage", pch=16)
legend("topright", legend=levels(cut(mean_log_cov, 7)),
       fill=brewer.pal(7, "Blues"))

# The trend being diagonal suggests the meth/GC content influences coverage
# This makes sense and is expected

# Lots of overplotting here, so dangerous to interpret
# But it looks like the motifs over-represented in WGBS are low GC content.
diff_cuts = cut(diffs, breaks = c(0.95, 0.98, 0.99, 1.01, 1.02, 1.09))
table(diff_cuts)
diff_cols = paste(brewer.pal(5, "Spectral")[diff_cuts], "80", sep="")
plot(seqs_mds, cex=0.5, col=diff_cols, main="Log2 Coverage")
legend("topright", legend=levels(diff_cuts),
       fill=brewer.pal(5, "Spectral"))


#####  Motif clustering using string distances  #####

# Distance matrix of motif strings
seqs_dist = stringdistmatrix(rownames(tt), method="hamming", useNames = TRUE)
# Hierarchical clustering
seqs_clust = hclust(seqs_dist)
# Let's plot to examine the dendrogram
plot(seqs_clust)
# On node 5 there's 4 families of seqeunces, let's group them
cuts = cutree(seqs_clust, k = NULL, h = 5)
table(cuts)
DNAStringSet(rownames(tt)[cuts == 1])

# Dom, these cut groups could be handed to your motif analysis (once we've taken care of the revcomp thing)

# The same hamming distance matrix can be used for MDS
seqs_mds = cmdscale(seqs_dist)
plot(seqs_mds, type="n")
text(seqs_mds, labels = rownames(seqs_mds), cex=0.6, col="grey50")

# Let's colour by GC content
gc_content = factor(str_count(rownames(tt), "G|C"))
gc_cols = brewer.pal(length(levels(gc_content)), "Spectral")[gc_content]
plot(seqs_mds, type="n")
text(seqs_mds, labels = rownames(seqs_mds), cex=0.6, col=gc_cols)


#####  Heatmaps  #####

# We can cluster the data using the coverage data, or enforce row order by string similarity

# e.g. clustering using a dist matrix and HClust on the coverage data of the most diff motifs
heatmap(log_motifs$all_enrich[rownames(tt), ])

motifs$all_enrich["TGACGTGT", ]

# Now we enforce the row order to be via the stingdistances dendrogram
# So now the HClust sorting is driven by motif hamming distances, not differences in coverage
# However, the heatmap colour is differences in coverage
seqs_clust = hclust(seqs_dist)
seqs_dend = as.dendrogram(seqs_clust)
heatmap(log_motifs$all_enrich[rownames(tt), ], Rowv = seqs_dend)


