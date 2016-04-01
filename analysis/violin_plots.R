doplot <- function(countdata, title) {
  bks = sapply(X=1:7, FUN=function(x) 10^x)
  
  ggplot(countdata) + 
    geom_violin(aes(x=project, y=num_mutations, fill=colour_project)) + 
    scale_y_log10(breaks=bks, labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    theme_grey() + annotation_logticks(sides="lr") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    xlab("Tumour Type") + ylab("Number of Mutations") + ggtitle(title) + 
    guides(fill=FALSE)
}

violin_plots <- function(indir='~/Desktop/hn.hpc/mounts/users/jdursi/interim_consensus/',
                         outdir='~/Desktop/',
                         snvfile='snv_mnv.summary.tsv',
                         indelfile='indel_normed.summary.tsv') {
  
  snv_counts <- read.table(paste0(indir,snvfile)); names(snv_counts) <- c("case","project","variant","num_mutations")
  indel_counts <- read.table(paste0(indir,indelfile)); names(indel_counts) <- c("case","project","variant","num_mutations")
  
  snv_ordered_tumour_types <- names(sort(tapply(snv_counts$num_mutations, snv_counts$project, FUN=median)))
  indel_ordered_tumour_types <- names(sort(tapply(indel_counts$num_mutations, indel_counts$project, FUN=median)))
  
  snv_counts$project <- factor(snv_counts$project, levels=snv_ordered_tumour_types)
  indel_counts$project <- factor(indel_counts$project, levels=indel_ordered_tumour_types)
  snv_counts$colour_project <- factor(snv_counts$project, levels=snv_ordered_tumour_types)
  indel_counts$colour_project <- factor(indel_counts$project, levels=snv_ordered_tumour_types)

  doplot(snv_counts, "Interim Merge Set: Number of passed SNVs")
  ggsave(paste0(outdir,"interim_merge_snv_violins.pdf"), height=10, width=12)
  
  doplot(indel_counts, "Interim Merge Set: Number of passed Indels")
  ggsave(paste0(outdir,"interim_merge_indel_violins.pdf"), height=10, width=12)
}