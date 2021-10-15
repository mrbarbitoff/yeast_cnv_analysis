setwd("/media/array/yeast_proj/plasmid_mutants")

library(ggplot2)
library(reshape2)
library(cowplot)
library(ggsci)
library(RColorBrewer)
library(ggrepel)


samples <- read.table('sample_map.csv', sep='\t', header=T)
head(samples)
ordered <- samples[c(1:9, 100, 10:99), ]
ordered_2 <- samples[c(1:10, 100, 11:99), ]

chrII <- read.table('qmap_II/binned_coverage.tsv', sep='\t', header=T, comment.char = '')
chrIIm <- melt(chrII, id.vars='X.Position..bp.')
head(chrIIm)

chrIIm$gene = rep(ordered$Gene, each=400)
chrIIm$mutant = rep(ordered$Mutant, each=400)

chrIIa <- aggregate(value~X.Position..bp.+gene+mutant, chrIIm, mean)

plot_II <- ggplot(chrIIa, aes(x=X.Position..bp., y=value, col=mutant)) +
  geom_line(lwd=1) + theme_bw() + 
  theme(panel.grid=element_blank()) + 
  facet_wrap(~gene)
print(plot_II)


chrIV <- read.table('qmap_IV/binned_coverage.tsv', sep='\t', header=T, comment.char = '')
chrIVm <- melt(chrIV, id.vars='X.Position..bp.')
head(chrIVm)

chrIVm$gene = rep(ordered$Gene, each=400)
chrIVm$mutant = rep(ordered$Mutant, each=400)

chrIVa <- aggregate(value~X.Position..bp.+gene+mutant, chrIVm, mean)

plot_IV <- ggplot(chrIVa, aes(x=X.Position..bp., y=value, col=mutant)) +
  geom_line(lwd=1) + theme_bw() + 
  theme(panel.grid=element_blank()) + 
  facet_wrap(~gene)
print(plot_IV)

plot_grid(plot_II, plot_IV, nrow=2)

total_agg <- rbind(chrIIa, chrIVa)
total_agg$chrom <- c(rep('II', nrow(chrIIa)), rep('IV', nrow(chrIVa)))

ggplot(total_agg, aes(x=X.Position..bp., y=value, col=mutant)) +
  geom_line(lwd=1) + theme_bw() + 
  theme(panel.grid=element_blank(), legend.position = 'bottom') + 
  facet_grid(vars(gene), vars(chrom), scales='free') +
  scale_color_jama() + xlab('Chromosomal coordinate') + ylab('Mean coverage')


### Panel c - WGS-based copies

wgs_35 <- read.table('pRSU1_WGS_copies.tsv', sep='\t', header=F)
wgs_45 <- read.table('pRS315_WGS_copies.tsv', sep='\t', header=F)

ordered_2$prsu <- wgs_35$V2
ordered_2$prs316 <- wgs_45$V2
ordered_2$copies = ifelse(ordered_2$Gene == 'sup35', ordered_2$prsu,
                          ordered_2$prs316)

wgs_data <- ordered_2[, c('SN', 'Gene', 'Mutant', 'Variant', 'Loss', 'copies')]
wgs_data$VariantLoss <- paste(wgs_data$Variant, wgs_data$Loss, sep='_')
wgs_data$mut_id <- sapply(strsplit(wgs_data$Variant, '-'), function(x) x[2])

f2c <- ggplot(wgs_data, aes(x=mut_id, y=copies, fill=Gene)) +
  geom_boxplot() + theme_bw() + guides(fill=F) +
  facet_wrap(~Gene, scales='free_x') +
  theme(panel.grid=element_blank()) +
  scale_fill_brewer(palette = "Accent") +
  scale_y_continuous(breaks = c(0, 3, 6, 9))
print(f2c)

ggplot(wgs_data, aes(x=VariantLoss, y=copies, fill=Gene)) +
  geom_boxplot() + theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1))

qpcr <- read.table('qPCR_plasmid_short.csv', sep=',', header=T)
f2d <- ggplot(qpcr, aes(x=mut_id, y=copies, fill=Gene)) +
  geom_boxplot() + theme_bw() + guides(fill=F) +
  facet_wrap(~Gene, scales='free_x') +
  theme(panel.grid=element_blank()) +
  scale_fill_brewer(palette = "Accent") + 
  scale_y_continuous(breaks=c(0, 3, 6, 9))
print(f2d)

f2cd <- plot_grid(f2c, f2d, nrow=2)
print(f2cd)

wgs_by_loss <- wgs_data[wgs_data$Gene == 'sup45', ]
head(wgs_by_loss)

ggplot(wgs_by_loss, aes(x=mut_id, y=copies, fill=Loss)) + 
  geom_boxplot() + theme_bw() + 
  theme(panel.grid=element_blank(), legend.position='bottom') +
  scale_fill_simpsons() + xlab('Allele') + ylab('Number of plasmid copies')


copy_agg <- aggregate(copies~mut_id+Gene, wgs_data, median)
copy_agg$qpcr <- aggregate(copies~mut_id+Gene, qpcr, median)$copies

f2e <- ggplot(copy_agg, aes(x=copies, y=qpcr)) + 
  geom_point(aes(fill=Gene), size=3, pch=21) + theme_bw() +
  scale_fill_brewer(palette = "Accent") +
  geom_smooth(method='lm') + guides(fill=F) +
  geom_text_repel(aes(label=mut_id))

plot_grid(f2cd, f2e, ncol=2)

cor.test(copy_agg$copies, copy_agg$qpcr)



### Figure 3

expr_data <- read.table('expression.csv', sep=',', header=T)
head(expr_data)
expr_data = expr_data[expr_data$Reference == 'ADH1', ]

ggplot(expr_data, aes(x=mut_id, y=Expression, fill=Gene)) +
  geom_boxplot() + theme_bw() + guides(fill=F) +
  theme(panel.grid=element_blank()) +
  facet_grid(vars(Target), vars(Gene), scales='free_x') +
  scale_fill_brewer(palette = "Accent") + 
  scale_y_continuous(limits=c(0, 9), breaks=c(0, 2, 4, 6, 8))




### Figure 4

chrom_qpcr <- read.table('/media/array/yeast_proj/chromosomal_mutatnts/1b_d1606_qPCR.csv',
                         sep='\t', header=T)
head(chrom_qpcr)

ggplot(chrom_qpcr, aes(x=mut_id, y=copies, fill=Gene)) + 
  geom_boxplot() + theme_bw() + guides(fill=F) +
  theme(panel.grid=element_blank()) + 
  facet_wrap(~Gene, scales='free_x', nrow=1) +
  scale_fill_brewer(palette = "Accent")


# 4B - coverage profiles
coverage_data <- read.table('/media/array/yeast_proj/chromosomal_mutatnts/coverage/all.bin.tsv', sep='\t', header=F)
colnames(coverage_data) = c('sample', 'contig', 'position', 'mean_cov', 'sd_cov')

cov_avg <- aggregate(mean_cov~sample, coverage_data, median)
coverage_data$norm_cov = coverage_data$mean_cov/rep(cov_avg$mean_cov, each=265)
coverage_data$norm_sd = coverage_data$sd_cov/rep(cov_avg$mean_cov, each=265)

coverage_data$lower = sapply(coverage_data$norm_cov - coverage_data$norm_sd,
                             function(x) max(x, 0.1))
coverage_data$upper = sapply(coverage_data$norm_cov + coverage_data$norm_sd,
                             function(x) min(x, 2.45))
coverage_data$norm_cov[coverage_data$norm_cov > 2.5] = 2.45
coverage_data$Gene = ifelse(coverage_data$sample %in% c('101', '102', '102-2',
                                                        '103', '104', '105',
                                                        '107'), 'sup45',
                            ifelse(coverage_data$sample == 'WT', 'wt', 'sup35'))

chr_bnd <- aggregate(position~contig, coverage_data, max)

ggplot(coverage_data, aes(x=position, y=norm_cov, col=Gene, fill=Gene)) + 
  geom_line() + geom_vline(data=chr_bnd, aes(xintercept=position), lwd=0.4, lty=2) +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.6, col=NA) +
  theme_bw() + facet_wrap(~sample, ncol=1) + 
  scale_y_continuous(limits=c(0, 2.5)) + 
  scale_color_brewer(palette = "Accent") + 
  scale_fill_brewer(palette = "Accent")

