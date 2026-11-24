#!/usr/bin/env Rscript
# Absolute native-network agreement for the seven Figure 3B metrics.
# Run from the project root. Existing main figures are not rebuilt or modified.
suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(grid)})
source('figures/shared/theme_paper.R')
source('figures/shared/manifest.R')
source('figures/shared/pairing_styles.R')
source('figures/shared/number_labels.R')

OUT <- 'figures/supplementary/figS_reuse_native_agreement'
LEVELS <- series_levels()
REUSE <- setNames(vapply(LEVELS, series_reuse, numeric(1)), LEVELS)
# Exact Figure 3 colour and shape mapping, including its reuse-specific fills.
RULE_HUE <- c(inference = '#D11F1F', random = '#2E7EBB')
FILLS <- c(inference.bijection = '#66000C', inference.cap14 = '#C0151A',
           inference.uncapped = '#FB6B4B', random.bijection = '#0F5AA3',
           random.cap14 = '#4292C6', random.uncapped = '#99C7E0')
SHAPES <- c(inference = 24, random = 21)
AML_RANGE <- sort(c(mf_num('aml.Diagnosis.achieved_reuse'),
                    mf_num('aml.Relapse.achieved_reuse')))
METRICS <- data.frame(
  stem = c('L1.rho','L1.topn_jaccard','L2.topn_jaccard','L3.genes',
           'L3.regions','L4.jaccard','L5.jaccard'),
  label = c('region-gene correlation\nSpearman rho',
            'top region-gene pairs\nJaccard index',
            'top TF-gene edges\nJaccard index',
            'eRegulon target genes\nmedian Jaccard index',
            'eRegulon target regions\nmedian Jaccard index',
            'signed eRegulon identities\nJaccard index',
            'TF identities\nJaccard index'),
  stringsAsFactors = FALSE)
D <- do.call(rbind, lapply(seq_len(nrow(METRICS)), function(j)
  do.call(rbind, lapply(LEVELS, function(l)
    do.call(rbind, lapply(c('inference','random'), function(r) {
      arm <- unname(series_arm(l, if (r == 'inference') 'link' else 'rand'))
      key <- paste('bmmc', METRICS$stem[j], arm, sep = '.')
      e <- mf_entry(key)
      data.frame(panel = letters[j], metric = METRICS$stem[j], rule = r,
                 level = l, reuse = REUSE[[l]], arm = arm,
                 agreement = mf_num(key), fill_key = paste(r,l,sep='.'),
                 manifest_key = key, source_file = e$source_file,
                 source_key = e$source_key, stringsAsFactors = FALSE)
    }))))))
stopifnot(nrow(D) == 42L, all(is.finite(D$agreement)),
          all(D$agreement >= 0 & D$agreement <= 1))
# Independently resolve every displayed value to its primary report.
for (i in seq_len(nrow(D))) {
  obj <- jsonlite::fromJSON(D$source_file[i], simplifyVector = FALSE)
  for (part in strsplit(D$source_key[i], '.', fixed = TRUE)[[1]]) obj <- obj[[part]]
  stopifnot(isTRUE(all.equal(as.numeric(obj), D$agreement[i], tolerance = 1e-12)))
}
trend <- vapply(split(D, interaction(D$metric,D$rule)), function(d)
  all(diff(d$agreement[order(d$reuse)]) < 0), logical(1))
message(sprintf('Primary-source check: 42/42 values match. Declining series: %d/14.',sum(trend)))
write.csv(D, paste0(OUT,'.csv'), row.names = FALSE)

make_panel <- function(j) {
  d <- D[D$metric == METRICS$stem[j],]
  p <- ggplot(d, aes(reuse,agreement,group=rule)) +
    annotate('rect', xmin=AML_RANGE[1], xmax=AML_RANGE[2],
             ymin=-Inf,ymax=Inf,fill=PAL[['context']],alpha=.16) +
    geom_hline(yintercept=1,colour='#111111',linewidth=.35) +
    geom_line(aes(colour=rule),linewidth=.4) +
    geom_point(aes(fill=fill_key,shape=rule),colour='white',size=2,stroke=.3) +
    scale_colour_manual(values=RULE_HUE,guide='none') +
    scale_fill_manual(values=FILLS,guide='none') +
    scale_shape_manual(values=SHAPES,guide='none') +
    scale_x_continuous(breaks=REUSE,labels=sprintf('%.3fx',REUSE),
                       expand=expansion(mult=c(.10,.10))) +
    labs(x='reuse',y=METRICS$label[j],tag=letters[j]) +
    theme_paper() +
    theme(axis.title=element_text(size=7),axis.text=element_text(size=6.3),
          plot.margin=margin(12,7,6,4),
          plot.tag=element_text(size=9,face='bold',hjust=0,vjust=1),
          plot.tag.position=c(0,1.03))
  if (j == 1) {
    p + scale_y_continuous(limits=c(.970,1.001),breaks=seq(.97,1,.01),
                           labels=lab_nolead,expand=expansion(mult=0))
  } else {
    p + scale_y_continuous(limits=c(0,1.02),breaks=seq(0,1,.25),
                           labels=lab_nolead,expand=expansion(mult=0))
  }
}
P <- lapply(seq_len(7),make_panel)
# The eighth grid position holds the common key; it is not a data panel.
K <- ggplot() +
  annotate('segment',x=.08,xend=.22,y=.80,yend=.80,colour=RULE_HUE[['inference']],linewidth=.4) +
  annotate('point',x=.15,y=.80,shape=24,fill=RULE_HUE[['inference']],colour='white',size=2.4,stroke=.3) +
  annotate('text',x=.27,y=.80,label='inference-assigned',hjust=0,size=2.6) +
  annotate('segment',x=.08,xend=.22,y=.62,yend=.62,colour=RULE_HUE[['random']],linewidth=.4) +
  annotate('point',x=.15,y=.62,shape=21,fill=RULE_HUE[['random']],colour='white',size=2.4,stroke=.3) +
  annotate('text',x=.27,y=.62,label='random within stratum',hjust=0,size=2.6) +
  annotate('rect',xmin=.08,xmax=.22,ymin=.39,ymax=.49,fill=PAL[['context']],alpha=.16) +
  annotate('text',x=.27,y=.44,label='pAML reuse range',hjust=0,size=2.6) +
  annotate('segment',x=.08,xend=.22,y=.26,yend=.26,colour='#111111',linewidth=.35) +
  annotate('text',x=.27,y=.26,label='perfect agreement',hjust=0,size=2.6) +
  coord_cartesian(xlim=c(0,1),ylim=c(0,1),expand=FALSE) + theme_void()
FIG <- wrap_plots(c(P,list(K)),ncol=2)
save_paper(FIG,OUT,w=COL2,h=200*MM,dpi=600)
mf_write_keys(OUT)

CAPTION <- paste0(
 'Supplementary Figure Sx. Agreement with the native network across RNA reuse levels.\n\n',
 'BMMC networks inferred under inferred and random pairing were compared with the ',
 'SCENIC+ network inferred from native same-cell RNA and ATAC measurements at reuse ratios of ',
 paste(sprintf('%.3f',REUSE),collapse=', '),'. The eligible RNA population, ATAC input and ',
 'within-donor-cell-type multiplicity distributions were matched between pairing rules at each reuse level. ',
 '(a) Spearman correlation of region-gene correlation coefficients across candidate pairs. ',
 '(b) Jaccard index of the top 100,000 region-gene pairs ranked by importance. ',
 '(c) Jaccard index of the top 100,000 TF-gene edges ranked by importance. ',
 '(d,e) Median target-gene and target-region Jaccard indices across shared eRegulons. ',
 '(f,g) Jaccard indices of signed eRegulon identities and TF identities. ',
 'Red triangles indicate inferred pairing; blue circles indicate random pairing within strata. ',
 'Marker fills become lighter at higher reuse, following the main Figure 3 palette. ',
 'The grey band marks the reuse range of the corrected pAML construction (',
 sprintf('%.4f',AML_RANGE[1]),' at diagnosis and ',sprintf('%.4f',AML_RANGE[2]),
 ' at relapse). Black horizontal lines indicate perfect agreement. ',
 'Panel a uses an expanded vertical scale from 0.97 to 1.00; panels b-g share a 0-1 scale. ',
 'Agreement declined across the three tested reuse levels for both rules on all seven metrics. ',
 'Each point represents one construction and network inference result; no repeated-construction ',
 'uncertainty intervals are shown. Connecting lines are visual guides between tested constructions. ',
 'The native network is a computational reference, not experimentally established regulatory ground truth.\n')
writeLines(CAPTION,paste0(OUT,'_caption.txt'),useBytes=TRUE)
