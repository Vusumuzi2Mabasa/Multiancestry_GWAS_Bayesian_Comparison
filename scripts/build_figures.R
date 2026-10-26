# Rebuild publication figures and derived tables from saved simulation outputs.
# No simulation is rerun and no package is installed by this script.
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly=FALSE)
script <- sub('^--file=','',args[grepl('^--file=',args)])
script <- gsub('~+~',' ',script,fixed=TRUE)
root <- normalizePath(file.path(dirname(script),'..'),mustWork=TRUE)
out <- file.path(root,'report');dir.create(out,showWarnings=FALSE)
figdir <- file.path(root,'results','figures');dir.create(figdir,showWarnings=FALSE)
stopifnot(trimws(readLines(file.path(root,'results','status.txt')))=='COMPLETE_WITH_MANC')
d <- read.delim(file.path(root,'results','replicate_metrics.tsv'),check.names=FALSE)
s <- read.delim(file.path(root,'results','summary.tsv'),check.names=FALSE)
pd <- read.delim(file.path(root,'results','paired_differences.tsv'),check.names=FALSE)
v <- read.delim(file.path(root,'results','variant_results.tsv'),check.names=FALSE)
cfgfile <- file.path(root,'results','config.json')
if(!file.exists(cfgfile))cfgfile<-file.path(root,'results','study_config.json')
cfg <- jsonlite::fromJSON(cfgfile)
stopifnot(cfg$seed==20261026,cfg$replicates==20,cfg$n_eur==2000,
          identical(as.integer(cfg$n_afr),c(250L,500L,1000L,2000L)))
prov <- jsonlite::fromJSON(file.path(root,'tools','provenance.json'))
stopifnot(nrow(d)==1200,nrow(v)==3840,all(table(d$method)==240))
primary <- c('Manc_COJO','Bayes_rho0_PIP0.95','Bayes_rho0.8_PIP0.95')
labels <- c(Manc_COJO='Manc COJO',Bayes_rho0_PIP0.95='Bayesian rho 0',
            Bayes_rho0.8_PIP0.95='Bayesian rho 0.8')
colors <- c('Manc COJO'='#286344','Bayesian rho 0'='#246b9e','Bayesian rho 0.8'='#b15c22')
z <- s[s$method%in%primary,]
z$Model <- factor(unname(labels[z$method]),levels=unname(labels))
z$Scenario <- factor(z$scenario,levels=c('shared','heterogeneous','null'),
                     labels=c('Shared effects','Different effects','Null control'))
theme_set(theme_bw(base_size=12,base_family='Arial')+
            theme(panel.grid.minor=element_blank(),legend.position='bottom',
                  legend.title=element_blank(),strip.background=element_rect(fill='#edf2f5'),
                  plot.title=element_text(face='bold',color='black',size=13),
                  plot.subtitle=element_text(size=10),plot.caption=element_text(hjust=0,size=9)))
saveplot <- function(name,p,h=3.7) {
  ggsave(file.path(figdir,name),p,width=6.75,height=h,units='in',dpi=300,bg='white')
}
zr <- z[z$scenario!='null',]
zr$lo<-pmax(0,zr$recall_mean-1.96*zr$recall_mcse)
zr$hi<-pmin(1,zr$recall_mean+1.96*zr$recall_mcse)
p <- ggplot(zr,aes(n_afr,recall_mean,color=Model,group=Model))+
  geom_line(linewidth=.7)+geom_errorbar(aes(ymin=lo,ymax=hi),width=45,linewidth=.4)+
  geom_point(size=2)+facet_wrap(~Scenario,nrow=1)+
  scale_color_manual(values=colors)+scale_x_continuous(breaks=cfg$n_afr)+
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.25))+
  labs(x='Smaller AFR labelled cohort size',y='Exact causal SNP recall',
       title='Recovery of the planted variants',
       subtitle='20 replicates per point    Bayesian PIP at least 0.95    COJO p below 5e-8')
saveplot('01_signal_recovery.png',p)
zr$rmse<-sqrt(zr$mse_afr_mean)
zr$rmse_se<-zr$mse_afr_mcse/(2*pmax(zr$rmse,1e-12))
p <- ggplot(zr,aes(n_afr,rmse,color=Model,group=Model))+
  geom_line(linewidth=.7)+geom_errorbar(aes(ymin=pmax(0,rmse-1.96*rmse_se),ymax=rmse+1.96*rmse_se),width=45,linewidth=.4)+
  geom_point(size=2)+facet_wrap(~Scenario,nrow=1)+scale_color_manual(values=colors)+
  scale_x_continuous(breaks=cfg$n_afr)+labs(x='Smaller AFR labelled cohort size',y='AFR effect RMSE across eight SNPs',
  title='Effect estimation in the smaller cohort',subtitle='Lower values indicate smaller error    Error bars use an approximate delta method')
saveplot('02_effect_error.png',p)
za<-z[z$scenario=='heterogeneous',]
p <- ggplot(za,aes(n_afr,afr_only_recall_mean,color=Model,group=Model))+
  geom_line(linewidth=.7)+geom_point(size=2)+
  geom_errorbar(aes(ymin=pmax(0,afr_only_recall_mean-1.96*afr_only_recall_mcse),
    ymax=pmin(1,afr_only_recall_mean+1.96*afr_only_recall_mcse)),width=45,linewidth=.4)+
  scale_color_manual(values=colors)+scale_x_continuous(breaks=cfg$n_afr)+
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.25))+
  labs(x='Smaller AFR labelled cohort size',y='AFR only variant recovery',title='Recovery of rsSIM006',
       subtitle='True effect is zero in EUR and positive 0.28 in AFR')
saveplot('03_afr_specific_recovery.png',p,h=3.1)
vs<-v[v$scenario=='heterogeneous' & v$SNP=='rsSIM006',]
vs$Model<-ifelse(vs$method=='Bayes_rho0','Bayesian rho 0','Bayesian rho 0.8')
bias<-do.call(rbind,lapply(split(vs,interaction(vs$n_afr,vs$Model,drop=TRUE)),function(a)
 data.frame(n_afr=a$n_afr[1],Model=a$Model[1],mean=mean(a$beta_afr),
            mcse=sd(a$beta_afr)/sqrt(nrow(a)))))
write.table(bias,file.path(root,'results','afr_specific_posterior_means.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
p<-ggplot(bias,aes(n_afr,mean,color=Model,group=Model))+
 geom_hline(yintercept=.28,linetype=2,color='black')+geom_line(linewidth=.7)+geom_point(size=2)+
 geom_errorbar(aes(ymin=mean-1.96*mcse,ymax=mean+1.96*mcse),width=45,linewidth=.4)+
 scale_color_manual(values=colors)+scale_x_continuous(breaks=cfg$n_afr)+
 labs(x='Smaller AFR labelled cohort size',y='Mean posterior AFR effect',
      title='Shrinkage of the ancestry specific effect',subtitle='Dashed line is the true AFR effect    Bayesian models only')
saveplot('04_afr_effect_shrinkage.png',p,h=3.1)

# A separate summary of discovery thresholds, with counts rather than rounded zeros.
null <- d[d$scenario=='null',]
nullcounts <- aggregate(any_false_positive~n_afr+method,data=null,sum)
write.table(nullcounts,file.path(root,'results','null_event_counts.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
