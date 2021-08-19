#####  Dependencies  #####

library(GenomicRanges)
library(RColorBrewer)
#library(vioplot)


#####  Paths  #####

path_base = getwd()
# Assumes your ./data symlink points to hb-stopwatch/work/cpgberus
path_cpgberus = file.path(path_base, "data")
path_cpgberus_covs = file.path(path_cpgberus, "04_parse_bismark_covs")
path_cpgberus_epic = file.path(path_cpgberus, "04_parse_methylationEPIC")
path_out = file.path(path_cpgberus, "11_methepic_vs_emseq_wgbs")

#####  Load data  #####

load(file=file.path(path_cpgberus_covs, "grch38p13_combined_covs_epic-overlap_grl.RData"))

#####  Sample  #####

# Combine meth and unmeth coverage into total coverage
for(s in names(covs_epic)) {
        covs_epic[[s]]$cov = covs_epic[[s]]$meth_cov + covs_epic[[s]]$unmeth_cov
        print(s)
        print(summary(covs_epic[[s]]$cov))
}

print(sapply(covs_epic, function(s) (median(s$cov))))
# WR025V1E WR025V1W WR025V9E WR025V9W WR069V1E WR069V1W WR069V9E WR069V9W 
# 23       10       13        7       38       12       11        6

# Each library has difference median coverage, but WGBS always lower than partner EM-seq lib
# Foreach lib pair, get the >= 75% quantile coverage for WGBS, capture all of these CpG,
# and for EM-seq, sample the same number of CpG sites as WGBS over that min-coverage.

set.seed(42)

covs_sample = list()
prefixes = unique(substr(names(covs_epic), 1, 7))

for(prefix in prefixes) {
        s_bis = paste(prefix, "W", sep="")
        s_em = paste(prefix, "E", sep="")
        
        cov_min = quantile(covs_epic[[s_bis]]$cov, 0.75)
        covs_sample[[s_bis]] = covs_epic[[s_bis]][covs_epic[[s_bis]]$cov >= cov_min, ]
        # Now sample the same number of rows from EM-seq
        x = covs_epic[[s_em]][covs_epic[[s_em]]$cov >= cov_min, ]
        covs_sample[[s_em]] = sample(x, size=length(covs_sample[[s_bis]]))
}


#####  Plots  #####


jet.colors = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
corr_lms = list()

pdf(file.path(path_out, "NGS_EPIC_correlation_plots_Quantile75.pdf"))
for(s in names(covs_sample)) {
        #summary(covs_epic[[s]]$cov)
        #cov_factor = cut(covs_epic[[s]]$cov, breaks=c(1, 5, 10, 20, 25, 30, 50, max(covs_epic[[s]]$cov)))
        #cov_cols = brewer.pal(length(levels(cov_factor)), "Spectral")[cov_factor]
        corr_lms[[s]] = lm(covs_sample[[s]]$meth_pct / 100 ~ covs_sample[[s]]$EPIC)
        
        smoothScatter(x = covs_sample[[s]]$EPIC,
                      y = covs_sample[[s]]$meth_pct / 100, nbin = 256,
                      main=paste(s, ", median:", median(covs_sample[[s]]$cov),
                                 ", nCpG:", length(covs_sample[[s]]), sep=""),
                      xlab="EPIC (beta)",
                      ylab=paste(s, "(beta)"),
                      colramp=jet.colors)
        abline(a=0, b=1, col="white", lty=2, lwd=2)
        abline(corr_lms[[s]], lwd=2, col="grey60")
}
dev.off()

df = data.frame("Median coverage"=sapply(covs_sample, function(x) median(x$cov)),
                "75% percentile coverage"=sapply(covs_sample, function(x) quantile(x$cov, 0.75)),
                "n"=sapply(corr_lms, function(x) nrow(x$model)),
                "intercept"=sapply(corr_lms, function(x) coefficients(x)[1]),
                "coefficient"=sapply(corr_lms, function(x) coefficients(x)[2]),
                "residual sum of squares"=sapply(corr_lms, deviance))

write.csv(df, file=file.path(path_out, "Coverage_and_lm_stats.csv"))

# Now with intercept == 0

corr_lms = list()

pdf(file.path(path_out, "NGS_EPIC_correlation_plots_Quantile75_0intercept.pdf"))
for(s in names(covs_sample)) {
        #summary(covs_epic[[s]]$cov)
        #cov_factor = cut(covs_epic[[s]]$cov, breaks=c(1, 5, 10, 20, 25, 30, 50, max(covs_epic[[s]]$cov)))
        #cov_cols = brewer.pal(length(levels(cov_factor)), "Spectral")[cov_factor]
        corr_lms[[s]] = lm(covs_sample[[s]]$meth_pct / 100 ~ 0 + covs_sample[[s]]$EPIC)
        
        smoothScatter(x = covs_sample[[s]]$EPIC,
                      y = covs_sample[[s]]$meth_pct / 100, nbin = 256,
                      main=paste(s, ", median:", median(covs_sample[[s]]$cov),
                                 ", nCpG:", length(covs_sample[[s]]), sep=""),
                      xlab="EPIC (beta)",
                      ylab=paste(s, "(beta)"),
                      colramp=jet.colors)
        abline(a=0, b=1, col="white", lty=2, lwd=2)
        abline(corr_lms[[s]], lwd=2, col="grey60")
}
dev.off()

df = data.frame("Median coverage"=sapply(covs_sample, function(x) median(x$cov)),
                "75% percentile coverage"=sapply(covs_sample, function(x) quantile(x$cov, 0.75)),
                "n"=sapply(corr_lms, function(x) nrow(x$model)),
                "intercept"=0,
                "coefficient"=sapply(corr_lms, function(x) coefficients(x)[1]),
                "residual sum of squares"=sapply(corr_lms, deviance))

write.csv(df, file=file.path(path_out, "Coverage_and_lm_stats_0intercept.csv"))


pdf(file.path(path_out, "NGS_EPIC_correlation_plots_Quantile75_TypeI.pdf"))
for(s in names(covs_sample)) {
        
        x = covs_sample[[s]][covs_sample[[s]]$Type == "I", ]
        
        res = lm(x$meth_pct / 100 ~ x$EPIC)
        
        smoothScatter(x = x$EPIC,
                      y = x$meth_pct / 100, nbin = 256,
                      main=paste(s, ", TypeI only, median:", median(x$cov),
                                 ", nCpG:", length(x), sep=""),
                      xlab="EPIC (beta)",
                      ylab=paste(s, "(beta)"),
                      colramp=jet.colors)
        abline(a=0, b=1, col="white", lty=2, lwd=2)
        abline(res, lwd=2, col="grey60")
}
dev.off()

# Half-baked playing around with stuff
if(FALSE) {
        
        #s=names(covs_epic)[6]
        #vioplot(covs_epic[[s]][covs_epic[[s]]$IlmnID %in% pan_tissue_imprinted]$meth_pct / 100, ylim=c(0, 1))
        
        half_meth_ps = lapply(names(covs_epic), function(s) {
                is_half = covs_epic[[s]]$EPIC > 0.45 & covs_epic[[s]]$EPIC < 0.55
                covs_epic[[s]]$IlmnID[is_half]
        })
        half_meth_ps = Reduce(intersect, half_meth_ps)
        
        sample_names = paste(rep(prefixes, each=3), c("E", "W", "I"), sep="")
        
        half_meth_m = matrix(NA_real_,
                             nrow=length(half_meth_ps),
                             ncol=length(prefixes) * 3,
                             dimnames = list(half_meth_ps, sample_names))
        
        for(s in names(covs_epic)) {
                
                is_half = covs_epic[[s]]$IlmnID %in% half_meth_ps
                is_high_cov = covs_epic[[s]]$cov >=10
                is_this = is_half & is_high_cov
                
                half_meth_m[covs_epic[[s]]$IlmnID[is_this], s] = covs_epic[[s]]$meth_pct[is_this] / 100
                half_meth_m[covs_epic[[s]]$IlmnID[is_this], paste(substr(s, 1, 7), "I", sep="")] = covs_epic[[s]]$EPIC[is_this]
        }
        
        
        boxplot(half_meth_m, las=2, cex.axis=0.8)
        abline(h=0.5, col="darkorange")
        
        smoothScatter(x = covs_epic[[s]]$EPIC,
                      y = covs_epic[[s]]$meth_pct, nbin = 512,
                      main=paste(s, ", median:", median(covs_epic[[s]]$cov), ", 3Q:", quantile(covs_epic[[s]]$cov, 0.75), sep=""),
                      xlab="EPIC (beta)",
                      ylab=paste(s, "(beta)"),
                      colramp=jet.colors)
        
        
        
        
        covs_epic[[s]]$cov = covs_epic[[s]]$meth_cov + covs_epic[[s]]$unmeth_cov
        #summary(covs_epic[[s]]$cov)
        cov_factor = cut(covs_epic[[s]]$cov, breaks=c(1, 5, 10, 20, 25, 30, 50, max(covs_epic[[s]]$cov)))
        cov_cols = brewer.pal(length(levels(cov_factor)), "Spectral")[cov_factor]
        
        
        plot(covs_epic[[s]]$EPIC,
             covs_epic[[s]]$meth_pct, pch=".",
             col=cov_cols,
             main=paste(s, ", median:", median(covs_epic[[s]]$cov), ", 3Q:", quantile(covs_epic[[s]]$cov, 0.75), sep=""),
             xlab="EPIC (beta)",
             ylab=paste(s, "(beta)")
        )
        
}
