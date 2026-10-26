# Rebuild the detailed Word report from the saved simulation outputs.
# No simulation is rerun and no package is installed by this script.
suppressPackageStartupMessages({library(officer);library(flextable);library(ggplot2)})
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

# officer's bundled template lacks Word Title and Subtitle styles. Add them
# in R before opening the authoring document so the title uses a real style.
template<-tempfile(fileext='.docx');print(read_docx(),target=template)
td<-tempfile('word_template_');dir.create(td);utils::unzip(template,exdir=td)
sx<-xml2::read_xml(file.path(td,'word','styles.xml'))
ns0<-xml2::xml_ns(sx)
for(st in xml2::xml_find_all(sx,'.//w:style',ns0)) {
 id<-xml2::xml_attr(st,'styleId')
 if(id %in% c('Normal','Titre1','Titre2','Titre3')) {
  xml2::xml_remove(xml2::xml_find_all(st,'./w:rPr',ns0))
  size<-if(id=='Normal')'22' else if(id=='Titre1')'30' else '24'
  bold<-if(id=='Normal')'' else '<w:b/>'
  rp<-xml2::read_xml(paste0('<w:rPr xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"><w:rFonts w:ascii="Arial" w:hAnsi="Arial"/>',bold,'<w:color w:val="000000"/><w:sz w:val="',size,'"/></w:rPr>'))
  xml2::xml_add_child(st,rp)
  xml2::xml_remove(xml2::xml_find_all(st,'./w:pPr',ns0))
  pptext<-if(id=='Normal')'<w:pPr xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"><w:spacing w:after="120" w:line="264" w:lineRule="auto"/></w:pPr>' else paste0('<w:pPr xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"><w:keepNext/><w:keepLines/><w:spacing w:before="180" w:after="100"/><w:outlineLvl w:val="',as.integer(sub('Titre','',id))-1,'"/></w:pPr>')
  xml2::xml_add_sibling(xml2::xml_find_first(st,'./w:rPr',ns0),xml2::read_xml(pptext),.where='before')
 }
}
for(item in list(c('Title','52'),c('Subtitle','28'))) {
 node<-xml2::read_xml(paste0('<w:style xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main" w:type="paragraph" w:styleId="',item[1],'">',
 '<w:name w:val="',item[1],'"/><w:basedOn w:val="Normal"/><w:next w:val="Normal"/><w:qFormat/>',
 '<w:pPr><w:keepNext/><w:spacing w:before="0" w:after="200"/></w:pPr>',
 '<w:rPr><w:rFonts w:ascii="Arial" w:hAnsi="Arial"/><w:color w:val="000000"/><w:sz w:val="',item[2],'"/></w:rPr></w:style>'))
 xml2::xml_add_child(sx,node)
}
xml2::write_xml(sx,file.path(td,'word','styles.xml'))
packed<-tempfile(fileext='.docx');zip::zipr(packed,list.files(td,all.files=TRUE,no..=TRUE),root=td)
doc <- read_docx(packed)
doc <- set_doc_properties(doc,title='Multiancestry GWAS simulation and Bayesian comparison',
                         subject='Methods and results of a paired synthetic pilot',
                         creator='Vusumuzi Mabasa',description='R generated report with methods before results')
doc <- body_set_default_section(doc,prop_section(
 page_size=page_size(width=8.27,height=11.69),
 page_margins=page_mar(top=.67,bottom=.65,left=.72,right=.72),
 footer_default=block_list(fpar(run_word_field('PAGE'),fp_p=fp_par(text.align='right'))),
 footer_even=block_list(fpar(run_word_field('PAGE'),fp_p=fp_par(text.align='right'))),
 footer_first=block_list(fpar(run_word_field('PAGE'),fp_p=fp_par(text.align='right')))))
txt <- function(t,style='Normal') {doc <<- body_add_par(doc,t,style=style)}
h1 <- function(t) txt(t,'heading 1')
h2 <- function(t) txt(t,'heading 2')
page <- function(t) {doc <<- body_add_break(doc);h1(t)}
cap <- function(t) {doc <<- body_add_fpar(doc,fpar(ftext(t,fp_text(bold=TRUE,color='black',font.size=10,font.family='Arial'))))}
picture <- function(name,height=3.7) {
 doc <<- body_add_img(doc,file.path(figdir,name),width=6.75,height=height)
}
tab <- function(dat,widths=NULL) {
 ft<-flextable(dat)
 ft<-font(ft,fontname='Arial',part='all');ft<-fontsize(ft,size=9.5,part='all')
 ft<-color(ft,color='black',part='all');ft<-bg(ft,bg='#183e55',part='header')
 ft<-color(ft,color='white',part='header');ft<-bold(ft,part='header')
 ft<-bg(ft,i=seq(2,nrow(dat),by=2),bg='#f1f5f8',part='body')
 ft<-border_remove(ft);border<-fp_border(color='#D9D9D9',width=.5)
 ft<-border_outer(ft,border=border,part='all');ft<-border_inner_h(ft,border=border,part='all')
 ft<-border_inner_v(ft,border=border,part='all');ft<-padding(ft,padding=5,part='all')
 ft<-valign(ft,valign='center',part='all');ft<-align(ft,align='center',part='all')
 ft<-align(ft,j=1,align='left',part='all')
 if(!is.null(widths))ft<-width(ft,j=seq_along(widths),width=widths) else ft<-autofit(ft)
 ft<-set_table_properties(ft,layout='fixed',width=1)
 ft<-paginate(ft,init=TRUE,hdr_ftr=TRUE)
 doc <<- body_add_flextable(doc,ft);txt('')
}
fmt <- function(x,digits=3)ifelse(is.finite(x),formatC(x,digits=digits,format='f'),'Not applicable')
ms <- function(x,se)paste0(fmt(x),' (',fmt(se),')')
primarytable <- function(sc) {
 a<-z[z$scenario==sc,];a<-a[order(a$n_afr,a$Model),]
 data.frame('AFR N'=a$n_afr,'Model'=as.character(a$Model),
   'Exact recall MCSE'=ms(a$recall_mean,a$recall_mcse),
   'Block recall'=fmt(a$causal_block_recall_mean),
   'AFR RMSE'=fmt(sqrt(a$mse_afr_mean)),'Mean FDP'=fmt(a$fdp_mean),check.names=FALSE)
}
metric <- function(sc,n,method,column) s[s$scenario==sc&s$n_afr==n&s$method==method,column][1]
diffrow <- function(sc,n,method,m)pd[pd$scenario==sc&pd$n_afr==n&pd$method==method&pd$metric==m,][1,]

txt('Multiancestry GWAS simulation and Bayesian comparison','Title')
txt('Detailed methodology and results','Subtitle')
txt('Prepared for Vusumuzi Mabasa')
txt('26 October 2026')
h1('Purpose of the study')
txt('This report evaluates how a smaller ancestry cohort affects the recovery of known genetic association signals and the estimation of their effects. It compares the actual Manc COJO program with a custom Bayesian hierarchical variable selection model using the same simulated datasets. The central question is whether borrowing information across cohorts helps when effects are shared and whether that borrowing becomes harmful when an effect is specific to the smaller cohort.')
txt('The report presents the methodology first, followed by results, interpretation, limitations and recommended extensions. The study is a small quantitative trait pilot designed to explain statistical behaviour and establish a reproducible workflow. It is not a genome wide performance benchmark or a claim of a new validated method.')
h2('Scope and terminology')
txt('EUR and AFR are labels for two artificial cohorts. They do not imply that a single LD pattern represents European or African genetic diversity. No AWI Gen, clinical or other real participant data were analysed. The simulation varies only selected features of LD, sample size and effect heterogeneity; it does not reproduce demographic history or admixture.')
txt('Throughout the report, rho refers to the prior correlation between ancestry specific SNP effects in the Bayesian model. It is different from the haplotype persistence parameter used to generate LD. Lower RMSE means better effect estimation. Recall means the fraction of planted causal variants selected, whereas block recall is a broader measure of locating the correct simulated region.')

page('Methodology')
h2('Study design and paired comparisons')
txt('The larger cohort contains 2000 people in every analysis. The smaller cohort contains 250, 500, 1000 or 2000 people. There are 20 independent simulation replicates for each of three effect scenarios. Four sample size settings multiplied by three scenarios and 20 replicates gives 240 paired analysis datasets. All methods analyse each dataset, so their differences can be calculated within the same replicate.')
cap('Table 1 Fixed simulation settings')
tab(data.frame(Setting=c('Larger cohort size','Smaller cohort sizes','Variants','Genomic organisation','Effect scenarios','Replicates','Random seed','Allele frequency','Haplotype persistence','Phenotype noise'),
 Value=c('2000','250 500 1000 2000','8 biallelic SNPs','Two independent blocks of four SNPs','Shared different and null','20 per scenario and sample size','20261026','0.30 in both simulated cohorts','EUR 0.80 and AFR 0.35','Independent normal noise with variance 1')),
 widths=c(2.1,4.65))
txt('Within one replicate, the EUR genotype and phenotype data are reused across the four smaller cohort sizes. The AFR datasets are nested subsets of the same simulated 2000 people. This reduces irrelevant sampling variation when studying increasing AFR sample size. Replicates are independent; sample size settings within a replicate are not independent.')
txt('These settings were fixed before the pilot results were inspected. The number of replicates is adequate for an initial workflow demonstration but insufficient for precise estimates of rare false positive events or a definitive comparison of statistical power.')

page('Genotype and phenotype generation')
h2('Generating different LD patterns')
txt('Each person receives two independently generated haplotypes. At the first SNP of each four SNP block, the A allele is sampled from a Bernoulli distribution with probability 0.30. At each subsequent SNP, the previous allele is retained with probability 0.80 in the EUR labelled cohort or 0.35 in the AFR labelled cohort. Otherwise a new allele is drawn with probability 0.30. The two haplotypes are added to give A allele dosages of 0, 1 or 2.')
txt('This process keeps the marginal allele frequency equal in expectation while producing stronger within block correlation in one cohort. The second block starts independently. This construction separates the consequences of LD and sample size from differences in minor allele frequency. It does not generate realistic chromosomes or attempt to represent all features of any ancestry.')
h2('Known genetic effects')
cap('Table 2 True per allele effects used in the simulation')
tab(data.frame(Scenario=c('Shared','Shared','Different','Different','Null'),SNP=c('rsSIM002','rsSIM006','rsSIM002','rsSIM006','All eight'),
 'EUR effect'=c('0.18','-0.16','0.18','0','0'),'AFR effect'=c('0.18','-0.16','0.09','0.28','0'),check.names=FALSE),
 widths=c(1.35,2.1,1.65,1.65))
txt('All SNPs not listed as nonzero in Table 2 have zero effects. In the different effect scenario, rsSIM006 is the deliberately AFR only signal, while rsSIM002 has a smaller effect in AFR. The study therefore tests both partial effect heterogeneity and a signal absent from the larger cohort.')
h2('Quantitative phenotype')
txt('The phenotype is the sum of centred genotype dosages multiplied by the true per allele effects, plus independent standard normal noise. Genotypes are centred before forming the genetic component. A phenotype is generated for the full AFR cohort before taking nested subsets, preserving the same individuals and their noise values across the sample size comparisons.')
txt('There are no missing genotypes, relatives, confounders, ancestry admixture or overlapping participants between cohorts. No adjustment covariates are needed in this deliberately homogeneous within cohort construction. Binary disease and survival outcomes are outside the model used here.')

page('Summary statistics and LD information')
h2('Marginal association estimates')
txt('Each SNP is tested separately within each cohort by ordinary least squares with an intercept. For a centred genotype column x and centred phenotype y, the marginal effect estimate is b = x transpose y divided by x transpose x. The standard error is computed from the residual variance of that marginal regression with n minus 2 degrees of freedom. Normal approximation two sided p values are supplied to Manc COJO.')
txt('Each summary statistics file contains SNP, A1, A2, freq, b, se, p and N in the order required by the software. A1 is always A, A2 is G and freq is the A allele frequency. All variants use the same synthetic chromosome and consistent positions, spaced 1000 bases apart. The effect scale is phenotype units per additional A allele, not a log odds ratio.')
h2('Reference genotypes and alignment')
txt('Each cohort has a PLINK BED BIM FAM reference dataset made from its own simulated genotypes. This gives in sample LD and removes external LD mismatch as a source of error. The Bayesian analysis receives the same LD information and exact genotype sums of squares. Both analyses therefore use the same underlying samples, alleles and correlation structure, although their internal statistical approximations differ.')
txt('BED encoding was validated by decoding the generated bytes and confirming that every genotype was recovered with the intended A allele dosage. The Manc COJO effect allele in selected output rows was also checked against the input A allele. This avoids comparing effects whose signs differ only because of allele coding.')
h2('Reconstructing the Bayesian likelihood')
txt('The Bayesian model does not use individual phenotypes. From marginal estimates and genotype sums of squares it reconstructs q = X transpose y. From the LD correlation matrix and genotype sums of squares it reconstructs Q = X transpose X. Conditional on the effect vector beta, q has a multivariate normal likelihood with mean Q times beta and covariance equal to residual variance times Q. The phenotype sum of squares is recovered from marginal regression estimates and standard errors. These identities were checked against direct calculations in a separate validation dataset.')
txt('Residual variance is estimated from the full eight SNP regression using the reconstructed sufficient statistics and n minus 9 degrees of freedom. That estimate is then held fixed in the Bayesian calculations. This is a plug in analysis: it does not integrate uncertainty in residual variance. The reconstructed linear model likelihood is appropriate for this small quantitative trait setting and cannot simply be transferred to logistic or survival GWAS summary statistics.')

page('Bayesian hierarchical variable selection')
h2('Variant inclusion and borrowing')
txt('Each SNP has a shared inclusion indicator gamma with prior probability 0.20. If gamma is zero, its effects are exactly zero in both cohorts. If gamma is one, the two effects follow a correlated normal prior. This explicit selection layer is what distinguishes the model from a dense normal shrinkage model in which every SNP retains a nonzero effect.')
txt('For an included SNP, the ancestry specific effect is a shared latent effect mu plus a cohort specific deviation delta. The variance of mu is rho times tau squared; each deviation has variance (1 - rho) multiplied by tau squared. Thus each ancestry effect has prior standard deviation tau and their prior correlation is rho. The effect scale is fixed at tau = 0.15 per allele.')
cap('Table 3 Bayesian prior settings and interpretation')
tab(data.frame(Parameter=c('Inclusion probability','Effect standard deviation','Primary correlation','Sensitivity correlation','Primary reporting threshold','Sensitivity reporting threshold'),
 Value=c('0.20 per SNP','0.15 per allele','rho 0.80','rho 0','PIP at least 0.95','PIP at least 0.50')),
 widths=c(3.35,3.4))
txt('The rho 0.80 model encourages similar effects in the two cohorts. The rho 0 sensitivity model removes effect size correlation but keeps the shared inclusion indicator. It is therefore not a fully separate ancestry baseline: evidence in one cohort can still influence whether a SNP is included for the pair of cohorts.')
h2('Posterior computation')
txt('All 256 possible inclusion configurations for eight SNPs are enumerated. For each configuration, Gaussian SNP effects are integrated analytically. The model evidence is combined with the configuration prior to obtain posterior model probabilities. Posterior inclusion probabilities and ancestry specific effect means are then averaged across configurations. No MCMC approximation is used.')
txt('The resulting posterior is exact conditional on the fixed prior settings and plug in residual variance. Hyperparameters are not learned from the data. The shared inclusion layer identifies the union of signals but cannot assign separate probabilities of an exactly zero effect within each ancestry. An AFR only effect can therefore be selected while its magnitude is pulled toward the larger cohort. The custom implementation is not MESuSiE or another established package.')

page('Manc COJO implementation')
h2('Actual software and selection settings')
txt('The comparison uses the official Manc COJO Windows executable from release v1.1.0. It was downloaded into the project tools directory, checked with its help command and run separately on all 240 datasets. A substitute stepwise implementation was not used. The source URL, byte count and SHA256 are recorded in tools provenance metadata.')
txt('For each run, the EUR and AFR PLINK reference prefixes are supplied in that order to the bfile option, and the corresponding summary statistics are supplied in the same order to the cojo file option. Stepwise selection uses the cojo slct mode, a p value threshold of 5 times 10 to the minus 8 and one computational thread. Other settings retain the executable defaults. Full commands, standard output, logs and selected SNP output files are saved locally for every dataset.')
txt('Manc COJO performs ancestry specific conditional and joint calculations and combines evidence in a meta analysis based selection procedure. Its output reports both cohort specific joint effects and a meta analysed joint effect. The analysis here uses cohort specific joint effects for effect estimation and the selected SNP list for recovery metrics. Unselected SNP effects are set to zero when defining a whole locus estimator.')
h2('Checks on completed runs')
txt('The audit confirmed that all 240 datasets have actual Manc COJO commands and captured output. Selected SNP effect alleles were A as expected. The joint meta analysis p values of selected SNPs were checked against the chosen selection threshold. Each method has exactly one result per scenario, sample size and replicate, and effect error measures are finite.')
h2('Different objectives and thresholds')
txt('Manc COJO is designed to identify approximately independent association signals. The Bayesian model expresses uncertainty over candidate variants. Selecting a correlated proxy can be a reasonable COJO result even if the proxy is not the exact planted causal SNP. Exact SNP recall is therefore reported together with a coarse block recall measure.')
txt('A PIP threshold of 0.95 and a p value threshold of 5 times 10 to the minus 8 are not equivalent. They do not guarantee the same false positive rate or false discovery rate. The current study does not calibrate thresholds to a common error rate. Any apparent difference in recovery must be interpreted with that limitation, rather than treated as a fair stand alone ranking of methods.')

page('Evaluation and Monte Carlo uncertainty')
h2('Recovery and false selections')
txt('Exact causal SNP recall is the number of selected planted causal SNPs divided by the number of planted causal SNPs. It is undefined in the null scenario. AFR only recall records whether rsSIM006 is selected in the different effect scenario. Block recall records the fraction of causal four SNP blocks containing at least one selected SNP. Because it accepts any selected variant in a block, block recall is deliberately permissive and is not a formal LD based signal matching rule.')
txt('False discovery proportion is the number of selected noncausal SNPs divided by the total selected SNPs, with zero assigned if nothing is selected. The report gives its average across replicates. A separate outcome records whether a replicate contains any false selection. Neither outcome should be confused with a per SNP genome wide type I error rate.')
h2('Effect estimation')
txt('For each cohort, mean squared error averages squared differences between estimated and true effects across all eight SNPs. RMSE is the square root of the mean MSE across replicates. Bayesian effects are posterior means for every SNP, regardless of the PIP threshold used to report discoveries. Manc COJO effects are joint estimates for selected SNPs and zero otherwise. This makes a whole locus comparison possible, but the methods have different selection and shrinkage behaviour.')
txt('MSE restricted to true causal SNPs is also saved in the numerical outputs. A low overall RMSE can arise from strong shrinkage at six null SNPs even if a true ancestry specific effect is biased. For that reason, the report examines the posterior mean at the planted AFR only variant separately.')
h2('Uncertainty and pairing')
txt('For a replicate level outcome, Monte Carlo standard error is the standard deviation across 20 replicates divided by the square root of 20. Paired differences subtract the Manc COJO outcome from the Bayesian outcome within each identical dataset before calculating the mean difference and its Monte Carlo standard error. SNPs are not treated as independent simulation replicates.')
txt('Plot error bars show mean plus or minus 1.96 Monte Carlo standard errors, truncated to zero and one for proportions. These are rough simulation uncertainty summaries, not Bayesian credible intervals or reliable small sample confidence bounds. RMSE bars use an approximate delta method. For zero events in 20 null replicates, a one sided 95 percent binomial upper bound is about 13.9 percent, so a zero observed error rate does not establish calibration.')

page('Reproducibility and report generation')
h2('Software and execution')
txt('Simulation and Bayesian calculations use Python 3.12.14 and NumPy 2.3.5 in the existing Windows runtime. Manc COJO runs locally on Windows. The initial summary tables and plots were created in base R 4.3.3 through WSL. This Word report is generated by the saved R script using officer, flextable and ggplot2, with numeric values read directly from the saved results. No analysis was submitted to an HPC cluster.')
txt('The main simulation script saves genotype reference files, marginal association statistics, truth and LD matrices for each dataset. A separate runner adds actual Manc COJO results to the existing datasets without regenerating the Bayesian analysis. The report can therefore be rebuilt without rerunning all simulations. Source code, compact result tables, figures and the Word document are suitable for the portfolio repository; generated genotype inputs and downloaded executables are excluded.')
h2('Validation and traceability')
txt('The statistical checks cover BED genotype encoding, reconstruction of sufficient statistics, a closed form posterior calculation under an orthogonal null model and invariance when the two ancestry inputs are swapped. The completed output audit covers 240 datasets, five method and threshold combinations and 3840 Bayesian SNP result rows. Configuration files and the executable checksum record the settings used.')
txt('The simulation uses one fixed random seed with independent streams indexed by scenario and replicate. Within replicate sample size comparisons are nested as described earlier. Result tables retain replicate level outcomes and Monte Carlo standard errors rather than only rounded summary plots. The posterior SD columns are standard deviations, not endpoints of a credible interval.')
h2('How to read the following results')
txt('The primary results use PIP at least 0.95 for both Bayesian prior settings and the fixed Manc COJO selection threshold. Sensitivity results at PIP at least 0.50 are presented separately. A larger recall indicates better exact recovery under the chosen decision rule. A smaller RMSE indicates better effect estimation under the whole locus estimator definition. Neither metric alone establishes overall superiority.')

page('Results for shared effects')
txt('Both cohorts have the same two causal effects in this scenario. At AFR N equal to 250, the rho 0.80 model recovers 67.5 percent of exact causal SNPs on average, compared with 60.0 percent for rho 0 and 45.0 percent for Manc COJO. With increasing AFR sample size, recovery improves for all three methods. At AFR N equal to 2000, the corresponding recalls are 100.0, 100.0 and 97.5 percent.')
cap('Table 4 Shared effect results at the primary decision thresholds')
tab(primarytable('shared'),widths=c(.5,1.45,1.55,1.05,1.05,1.15))
txt('Exact recall is shown as the mean with Monte Carlo standard error in parentheses. Block recall and FDP are means across replicates. RMSE is calculated from the mean MSE. These thresholds do not represent matched error rates.')
h2('Interpretation')
txt('When effects are shared, the prior encouraging similar ancestry effects can use the larger cohort to stabilise estimates in the smaller cohort. At AFR N equal to 250, AFR RMSE is 0.039 for rho 0.80, 0.050 for rho 0 and 0.077 for Manc COJO. The benefit is smaller when the AFR cohort reaches 2000 people, consistent with less need to borrow information.')
dr<-diffrow('shared',250,'Bayes_rho0.8_PIP0.95','recall')
txt(paste0('The paired difference in exact recall between rho 0.80 and Manc COJO at AFR N equal to 250 is ',fmt(dr$difference),' with Monte Carlo standard error ',fmt(dr$mcse),'. This describes the chosen settings in this pilot; the comparison still combines different decision rules and should not be interpreted as a calibrated power advantage.'))

page('Signal recovery across sample sizes')
cap('Figure 1 Exact causal SNP recovery in the two nonnull scenarios')
picture('01_signal_recovery.png')
txt('Points are averages over 20 paired replicates. Error bars are approximate Monte Carlo uncertainty summaries. The two panels use the same model definitions, but the shared scenario has identical effects across cohorts whereas the different effect scenario contains an AFR only effect. Exact recall concerns the planted SNPs rather than all statistically equivalent proxies.')
h2('What the pattern means')
txt('Increasing the smaller cohort size improves the information available for distinguishing correlated variants and estimating their effects. The effect is particularly important when a signal is not present in the larger cohort: the large EUR sample cannot directly supply the missing AFR association evidence. A method that borrows across cohorts can still benefit from shared signals, but assumptions about similarity matter.')
txt('The larger gap between the Bayesian models and Manc COJO in the different effect scenario is not proof that the Bayesian model is universally better. The stepwise meta analysis criterion can dilute a signal present in only one cohort, and the Bayesian threshold is not matched to the stringent COJO p value rule. Both the selection objective and the error rate operating point contribute to the displayed differences.')

page('Results when effects differ across cohorts')
txt('In this scenario, rsSIM002 has a smaller AFR effect than EUR effect and rsSIM006 is present only in AFR. At AFR N equal to 250, exact recall is 30.0 percent for rho 0.80, 27.5 percent for rho 0 and 20.0 percent for Manc COJO. These aggregate recalls combine two signals with different architectures and should be read alongside the AFR only analysis.')
cap('Table 5 Different effect results at the primary decision thresholds')
tab(primarytable('heterogeneous'),widths=c(.5,1.45,1.55,1.05,1.05,1.15))
txt('Mean FDP is 0.050 for rho 0.80 and 0.025 for rho 0 at AFR N equal to 250, compared with zero observed for Manc COJO at that setting. Some gains in recall therefore coexist with more false selections. Twenty replicates leave substantial Monte Carlo uncertainty.')
h2('Interpretation')
txt('At AFR N equal to 250, rho 0 has a lower AFR RMSE of 0.078 than the stronger borrowing model at 0.086. The direction is the reverse of the shared effect scenario. This is the central caution: encouraging similar effects can reduce variance when the assumption is appropriate and increase bias when it is not. Manc COJO has RMSE 0.105 under the whole locus estimator definition, reflecting the combination of nonselection, effect estimation and its different threshold.')

page('Effect estimation and the borrowing tradeoff')
cap('Figure 2 AFR effect estimation error across all eight SNPs')
picture('02_effect_error.png')
txt('RMSE is the square root of mean replicate MSE across all eight SNPs. Error bars use an approximate delta method and describe Monte Carlo uncertainty, not uncertainty for a single participant or SNP. The selected estimator for Manc COJO assigns zero effects to unselected variants; Bayesian posterior means retain model averaged shrinkage.')
h2('What lower error does and does not imply')
txt('The Bayesian models estimate nonzero posterior means even when PIP falls below the reporting threshold. Their estimation and discovery decisions are therefore different operations. Lower whole locus error can be valuable, but it does not necessarily mean the method discovered more independent signals or isolated the causal variant more confidently.')
txt('The comparison between rho 0 and rho 0.80 is informative because the likelihood, inclusion prior, effect scale and computation are otherwise held fixed. Strong borrowing lowers error in the shared scenario and raises it in the different effect scenario. This supports investigating the amount of borrowing rather than assuming that more pooling always benefits underrepresented cohorts.')
txt('A practical follow up would allow more flexible shared and ancestry specific effect configurations and evaluate sensitivity to effect prior scale. Those changes would need prespecified evaluation and should not be tuned on these same 20 replicates to make one method appear better.')

page('The ancestry specific variant')
cap('Figure 3 Selection of the AFR only variant rsSIM006')
picture('03_afr_specific_recovery.png',height=3.1)
txt('At AFR N equal to 250, the rho 0 model selects the AFR only variant in 10 percent of replicates at PIP at least 0.95, while rho 0.80 and Manc COJO select it in none. By AFR N equal to 1000, both Bayesian models select it in 90 percent of replicates, compared with 5 percent for Manc COJO under the chosen settings.')
cap('Figure 4 Bayesian posterior mean for the AFR only effect')
picture('04_afr_effect_shrinkage.png',height=3.1)
txt('The true AFR effect is 0.28 and the EUR effect is zero. The dashed line marks the AFR truth. Strong borrowing pulls the AFR estimate further toward the larger cohort when the smaller cohort has limited information. Error bars summarise variation in the mean across replicates; they are not posterior credible intervals.')

page('Sensitivity to the discovery threshold')
txt('The PIP threshold controls which posterior signals are reported as discoveries. Reducing it from 0.95 to 0.50 raises recall but can increase false selections. It does not change the posterior mean effect estimates. Consequently, RMSE is identical between the two reporting thresholds within a Bayesian prior setting.')
cap('Table 6 Discovery sensitivity at AFR N equal to 250')
a<-s[s$n_afr==250&s$scenario!='null',]
a<-a[order(a$scenario,a$method),]
methodlabel<-function(m)ifelse(m=='Manc_COJO','Manc COJO',
  ifelse(grepl('rho0.8',m,fixed=TRUE),'Bayes rho 0.8','Bayes rho 0'))
sens<-data.frame(Scenario=ifelse(a$scenario=='shared','Shared','Different'),Model=methodlabel(a$method),
 Rule=ifelse(a$method=='Manc_COJO','p below 5e-8',ifelse(grepl('PIP0.95',a$method,fixed=TRUE),'PIP 0.95','PIP 0.50')),
 Recall=fmt(a$recall_mean),'Mean FDP'=fmt(a$fdp_mean),'Any false'=fmt(a$any_false_positive_mean),check.names=FALSE)
tab(sens,widths=c(.85,1.35,1.25,1.0,1.1,1.2))
txt('Any false is the proportion of replicates with at least one false selection. FDP is zero in a replicate where nothing is selected. The table demonstrates why comparing recall at arbitrary thresholds is not sufficient to conclude which procedure has better calibrated power.')
h2('Implication for the next benchmark')
txt('A stronger benchmark would use independent calibration simulations to select operating points at comparable error rates, then evaluate recovery on held out simulated datasets. Exact causal variant recovery and LD equivalent signal recovery should both be reported. Prior scale sensitivity and flexible heterogeneity models should be evaluated without choosing settings after seeing which ones look best.')

page('Null controls and quality checks')
txt('The null scenario has no true effects in either cohort. It measures false selection behaviour in these eight SNP loci, not genome wide family wise error. The primary Bayesian threshold and Manc COJO produced no false selections across the 20 null replicates at each AFR sample size. Some false selections appear at the lower PIP threshold.')
cap('Table 7 Null replicates with any false selection')
nr<-data.frame('AFR N'=cfg$n_afr,check.names=FALSE)
for(m in c('Manc_COJO','Bayes_rho0_PIP0.95','Bayes_rho0.8_PIP0.95','Bayes_rho0_PIP0.5','Bayes_rho0.8_PIP0.5')) {
 col<-if(m=='Manc_COJO')'Manc COJO' else paste(if(grepl('rho0.8',m,fixed=TRUE))'rho 0.8' else 'rho 0',if(grepl('0.95',m,fixed=TRUE))'PIP 0.95' else 'PIP 0.50')
 nr[[col]]<-vapply(cfg$n_afr,function(n)paste0(sum(null$any_false_positive[null$n_afr==n&null$method==m]),' / 20'),character(1))
}
tab(nr,widths=c(.55,1.15,1.25,1.25,1.25,1.3))
txt('Zero observed events in 20 replicates has a one sided 95 percent upper bound of approximately 13.9 percent for the probability of a replicate containing a false selection. Monte Carlo standard errors can be zero when no events are observed; that numerical zero does not imply certainty about the underlying rate.')
h2('Interpretation of null effect error')
txt('When Manc COJO selects nothing, its all zero effect estimate has zero error under an all zero truth. Bayesian posterior means generally remain slightly nonzero, so their null RMSE is positive even when no SNP reaches the discovery threshold. That difference follows from the estimator definitions and does not by itself imply a calibration failure.')
h2('Completed output audit')
txt('All 240 paired datasets have five method and threshold result rows. The audit verified 3840 Bayesian SNP rows, finite MSE values, valid PIP ranges, actual Manc COJO command and standard output files, consistent effect alleles and selected joint p values. Statistical unit checks passed for BED encoding, sufficient statistics, an analytic null posterior and ancestry swap invariance.')

page('Discussion and practical meaning')
h2('What this pilot supports')
txt('The pilot supports a conditional conclusion: borrowing information can help a smaller ancestry cohort when effects are similar, but stronger assumptions of similarity can hurt an effect that differs or is absent in the larger cohort. The smaller cohort sample size remains important even when a large external cohort exists. More data from the smaller cohort improve information about its own association pattern and reduce reliance on prior assumptions.')
txt('The results do not establish that the custom Bayesian model is superior to Manc COJO. It operates at a different discovery threshold, uses a different selection objective and estimates effects through model averaging. Its favourable whole locus RMSE partly reflects shrinkage and the treatment of unselected variants. The comparison is useful for understanding these choices, not for announcing a winning method.')
h2('Implications for African ancestry studies')
txt('In a real study, smaller African ancestry sample sizes can make association and fine mapping uncertain. Borrowing across cohorts may improve precision, but it should not force all effects to resemble those in the largest cohort. LD references must also be appropriate to the samples being analysed. The artificial AFR label here is insufficient to establish performance in any real African population, and the same modelling principles apply whenever one cohort is smaller or genetically heterogeneous.')
h2('Limitations that change the interpretation')
txt('The study contains only eight variants and two simple blocks. It uses clean in sample LD, identical expected allele frequencies, fixed prior parameters and a quantitative phenotype. There is no reference mismatch, admixture, missingness, imputation uncertainty, sample overlap or confounding. The plug in residual variance is treated as known in the posterior calculation. Twenty replicates provide coarse simulation uncertainty and little information about rare false positives.')
txt('The Bayesian shared inclusion indicator does not separately represent whether an effect is exactly absent in one ancestry. Rho 0 still shares inclusion evidence. Block recall is a coarse regional measure. No posterior credible set coverage, predictive performance, genome wide scalability or formal error matched power comparison has been evaluated.')

page('Recommended next study')
h2('Prespecified extensions')
txt('The next phase should increase the number of independent replicates, particularly under the null, and introduce larger loci with several realistic LD patterns. Sample size imbalance should be varied while retaining paired inputs across methods. Simulations should include shared effects, ancestry specific effects, opposite effects and different allele frequencies rather than varying only one favourable architecture.')
txt('Separate reference panels should be simulated to test finite reference sample size and LD mismatch. The evaluation should distinguish exact causal SNP recovery, LD equivalent signal recovery, false selections, effect bias, effect RMSE and uncertainty calibration. Discovery thresholds should be calibrated on independent data and then assessed on held out replicates.')
txt('An established Bayesian multiancestry fine mapping method such as MESuSiE could provide an additional benchmark. The custom model could be extended to ancestry specific inclusion states, learned effect correlations or robust mixtures, but each added feature would require independent validation. A new hierarchy is not automatically a new methodological contribution.')
h2('Decision from the current pilot')
txt('The workflow is ready to support a larger, prespecified methodological experiment. It should be described as a reproducible synthetic pilot demonstrating the benefits and risks of cross cohort borrowing. Its immediate value is the transparent connection between simulation truth, analysis choices and observed outcomes.')

page('References and reproducibility files')
txt('1. Wang X and colleagues. Multiancestry conditional and joint analysis applied to GWAS summary statistics. Nature Genetics 2026. https://doi.org/10.1038/s41588-026-02799-1')
txt('2. Manc COJO official tutorial. Input format, multiple cohort ordering and joint effect output. https://light156.github.io/multi-ancestry-COJO-docs/tutorial/')
txt('3. Manc COJO official command options. https://light156.github.io/multi-ancestry-COJO-docs/options/')
txt('4. Manc COJO official source and release downloads. https://github.com/light156/multi-ancestry-COJO')
txt('5. MESuSiE enables scalable and powerful multi ancestry fine mapping of causal variants in genome wide association studies. Nature Genetics 2024. https://doi.org/10.1038/s41588-023-01604-7')
h2('Files supporting this report')
txt('scripts/simulate_compare.py contains the simulation and custom Bayesian model. scripts/run_manc_existing.py runs the actual Manc COJO executable on existing inputs. scripts/validate.py and scripts/audit_results.py implement the statistical and completed output checks. scripts/report.R creates compact summary tables and the initial overview. scripts/build_word_report.R creates this Word document and its ggplot2 figures.')
txt('results/replicate_metrics.tsv preserves paired replicate outcomes. results/variant_results.tsv contains Bayesian PIPs, ancestry specific posterior means, posterior SDs and truth. results/summary.tsv and results/paired_differences.tsv contain the aggregate statistics and Monte Carlo standard errors. Additional report derived tables record null event counts and the AFR only posterior means.')
txt('The official Manc COJO release used here is v1.1.0. Its SHA256 is recorded in tools/provenance.json. The executable itself is excluded from the GitHub repository; a user can obtain it from the official release and verify it against that provenance. Generated genotype inputs and per run logs remain local and can be reproduced from the saved seed and code.')

outfile<-file.path(out,'Multiancestry_GWAS_Methodology_and_Results.docx')
print(doc,target=outfile)
if(Sys.getenv('WORD_REPORT_KEEP_RAW')=='1')file.copy(outfile,file.path(out,'.raw_before_style.docx'),overwrite=TRUE)
capture.output(sessionInfo(),file=file.path(out,'R_report_sessionInfo.txt'))
jsonlite::write_json(list(R=R.version.string,packages=setNames(
 lapply(c('officer','flextable','ggplot2','jsonlite','xml2','zip'),function(p)as.character(packageVersion(p))),
 c('officer','flextable','ggplot2','jsonlite','xml2','zip'))),
 file.path(root,'results','word_report_software.json'),pretty=TRUE,auto_unbox=TRUE)
message('Created ',outfile)
