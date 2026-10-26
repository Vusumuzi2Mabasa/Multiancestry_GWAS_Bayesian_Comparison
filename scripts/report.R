args <- commandArgs(trailingOnly=FALSE)
script <- sub('^--file=', '', args[grepl('^--file=',args)])
script <- gsub('~+~',' ',script,fixed=TRUE)
root <- normalizePath(file.path(dirname(script),'..'),mustWork=TRUE)
d <- read.delim(file.path(root,'results','replicate_metrics.tsv'),check.names=FALSE)
out <- file.path(root,'results')
capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
keys <- c('scenario','n_afr','method')
measures <- c('recall','fdp','any_false_positive','afr_only_recall',
              'causal_block_recall','mse_afr','mse_eur','mse_afr_causal','seconds')
groups <- split(seq_len(nrow(d)), interaction(d[keys],drop=TRUE))
summary <- do.call(rbind,lapply(groups,function(i) {
  z <- d[i,,drop=FALSE]; base <- z[1,keys,drop=FALSE];base$replicates<-nrow(z)
  for (m in measures) {
    x<-z[[m]];x<-x[is.finite(x)]
    base[[paste0(m,'_mean')]]<-if(length(x))mean(x) else NA_real_
    base[[paste0(m,'_mcse')]]<-if(length(x)>1)sd(x)/sqrt(length(x)) else NA_real_
  }
  base
}))
rownames(summary)<-NULL
write.table(summary,file.path(out,'summary.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
methods<-unique(as.character(summary$method))
cols<-setNames(c('#174a7e','#57a0ce','#b65520','#e9ae60','#3f8b58')[seq_along(methods)],methods)
plots <- function() {
  layout(rbind(1:3,4:6,c(7,7,7)),heights=c(1,1,.22))
  par(mar=c(4,4,3,1),oma=c(0,0,3,0))
  for(s in c('shared','heterogeneous','null')) {
    z<-summary[summary$scenario==s,,drop=FALSE]
    metric<-if(s=='null')'any_false_positive_mean' else 'recall_mean'
    plot(NA,xlim=c(250,2000),ylim=c(0,1),xlab='Smaller AFR-labelled cohort N',
      ylab=if(s=='null')'Replicates with any false selection' else 'Exact causal-SNP recall',main=s)
    for(m in methods) {a<-z[z$method==m,,drop=FALSE];a<-a[order(a$n_afr),]
      lines(a$n_afr,a[[metric]],type='b',pch=16,col=cols[m],lwd=2)}
  }
  for(s in c('shared','heterogeneous','null')) {
    z<-summary[summary$scenario==s,,drop=FALSE]
    plot(NA,xlim=c(250,2000),ylim=range(c(0,sqrt(z$mse_afr_mean)),na.rm=TRUE),
      xlab='Smaller AFR-labelled cohort N',ylab='AFR effect RMSE, all 8 SNPs',main=s)
    for(m in methods) {a<-z[z$method==m,,drop=FALSE];a<-a[order(a$n_afr),]
      lines(a$n_afr,sqrt(a$mse_afr_mean),type='b',pch=16,col=cols[m],lwd=2)}
  }
  mtext('Synthetic small-locus pilot: thresholds differ; no matched-error claim',outer=TRUE,cex=1.1)
  par(mar=c(0,0,0,0))
  plot.new();legend('center',legend=methods,col=cols[methods],lwd=2,pch=16,bty='n',ncol=2,cex=.8)
}
png(file.path(out,'comparison.png'),width=2100,height=1400,res=170);plots();dev.off()
png(file.path(out,'legend.png'),width=1400,height=350,res=150)
par(mar=c(0,0,0,0));plot.new();legend('center',legend=methods,col=cols[methods],lwd=3,pch=16,bty='n',cex=.9)
dev.off()
# Paired differences use replicate IDs within a scenario/sample size.
baseline<-if('Manc_COJO'%in%methods)'Manc_COJO' else 'Bayes_rho0_PIP0.95'
paired<-list();counter<-0
for(m in setdiff(methods,baseline)) {
 a<-d[d$method==m,];b<-d[d$method==baseline,]
 z<-merge(a,b,by=c('scenario','n_afr','replicate'),suffixes=c('.a','.b'))
 for(s in unique(z$scenario))for(n in unique(z$n_afr)) {
  zz<-z[z$scenario==s & z$n_afr==n,]
  for(metric in c('recall','fdp','afr_only_recall','causal_block_recall','mse_afr')) {
   x<-zz[[paste0(metric,'.a')]]-zz[[paste0(metric,'.b')]];x<-x[is.finite(x)]
   if(length(x)) {counter<-counter+1;paired[[counter]]<-data.frame(scenario=s,n_afr=n,
     method=m,baseline=baseline,metric=metric,difference=mean(x),
     mcse=if(length(x)>1)sd(x)/sqrt(length(x)) else NA,n=length(x))}
  }
 }
}
write.table(do.call(rbind,paired),file.path(out,'paired_differences.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
cat('Wrote summary tables, paired comparisons, comparison.png and legend.png\n')
status<-readLines(file.path(out,'status.txt'))
sink(file.path(out,'RESULTS.md'))
cat('# Synthetic multiancestry pilot results\n\nStatus: **',status,'**\n\n',sep='')
if(!'Manc_COJO'%in%methods)cat('Manc-COJO has not been run. The tables currently compare Bayesian prior settings only.\n\n')
cat('## Design and interpretation\n\n',
    'Eight SNPs, 2,000 EUR-labelled individuals and nested AFR-labelled samples of 250, 500, 1,000 and 2,000. ',
    'Three scenarios; ',max(summary$replicates),' replicates per setting. ',
    'These are stylised labels, not real ancestry data. See README.md for the complete model and limitations.\n\n',sep='')
cat('## Primary PIP >= 0.95 results\n\n',
    '| Scenario | AFR N | Method | Exact recall | AFR-only recall | AFR RMSE | Mean FDP |\n',
    '|---|---:|---|---:|---:|---:|---:|\n',sep='')
z<-summary[grepl('PIP0.95$',summary$method)|summary$method=='Manc_COJO',]
z<-z[order(z$scenario,z$n_afr,z$method),]
fmt<-function(x)ifelse(is.finite(x),sprintf('%.3f',x),'NA')
for(i in seq_len(nrow(z)))cat('| ',z$scenario[i],' | ',z$n_afr[i],' | ',z$method[i],
  ' | ',fmt(z$recall_mean[i]),' | ',fmt(z$afr_only_recall_mean[i]),' | ',
  fmt(sqrt(z$mse_afr_mean[i])),' | ',fmt(z$fdp_mean[i]),' |\n',sep='')
cat('\n## Limits\n\n',
    '- Recall counts exact planted variants. Correlated proxies may recover a region without recovering the exact SNP; block recall is also tabulated.\n',
    '- PIP and p-value thresholds are not matched for error rate. Do not infer method superiority directly from a recall difference.\n',
    '- Zero false selections in 20 null replicates has a one-sided 95% upper bound of about ',
    sprintf('%.1f%%',100*(1-.05^(1/max(summary$replicates)))),
    ' on the per-replicate probability; it does not establish genome-wide calibration.\n',
    '- Shrinkage can improve total effect RMSE while biasing ancestry-specific effects. Inspect per-variant estimates, not only the aggregate score.\n',
    '- The rho=0 model retains shared variant inclusion; it is not a fully independent-ancestry baseline.\n',
    '- In-sample LD, fixed effect priors and a small quantitative locus favour clean calculations; external LD, binary/survival outcomes and real population structure remain untested.\n\n',
    'Tables: summary.tsv (means and Monte Carlo SEs), paired_differences.tsv, replicate_metrics.tsv and variant_results.tsv. ',
    'Figure: comparison.png. No real participant data were used.\n',sep='')
sink()
