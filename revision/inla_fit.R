# Independent R-INLA fits to the whitened summary likelihood. No CAR prior.
args <- commandArgs(trailingOnly=TRUE)
suppressPackageStartupMessages(library(INLA))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(jsonlite))
root <- normalizePath('.', winslash='/')
source(file.path(root,'revision/inla_mixture.R'))
gate <- fromJSON(file.path(root,'results/revision/baseline_gate.json'))
stopifnot(isTRUE(gate$passed),isTRUE(gate$baseline_completed_before_inla))
inp <- fromJSON(args[1],simplifyVector=FALSE)
output <- args[2]
rows <- list()
for (i in seq_along(inp)) {
  z <- inp[[i]]; y <- unlist(z$y)
  B <- matrix(unlist(z$B),nrow=16,ncol=16,byrow=TRUE)
  Am <- B %*% rbind(diag(8),diag(8)); Ad <- B
  stk <- inla.stack(data=list(y=y),A=list(Matrix(Am,sparse=TRUE),Matrix(Ad,sparse=TRUE)),effects=list(data.frame(mu=1:8),data.frame(delta=1:16)))
  hp <- if (z$mode=='fixed') list(prec=list(initial=log(1/.075^2),fixed=TRUE)) else list(prec=list(initial=log(1/.075^2),prior='pc.prec',param=c(.2,.05)))
  form <- y ~ -1 + f(mu,model='iid',n=8,hyper=list(prec=list(initial=log(1/.15^2),fixed=TRUE))) + f(delta,model='iid',n=16,hyper=hp)
  lc <- inla.make.lincombs(mu=rbind(diag(8),diag(8)),delta=diag(16))
  integration <- if (z$mode=='fixed') list(strategy='gaussian',int.strategy='grid') else list(strategy='gaussian',int.strategy='user',int.design=cbind(seq(-2,28,length.out=201),c(.5,rep(1,199),.5)/200))
  integration$npoints <- 101L
  integration$cutoff <- 1e-12
  integration$diff.logdens <- 20
  start <- proc.time()[3]
  fit <- inla(form,family='gaussian',data=inla.stack.data(stk),lincomb=lc,
    control.family=list(hyper=list(prec=list(initial=0,fixed=TRUE))),
    control.predictor=list(A=inla.stack.A(stk),compute=FALSE),
    control.inla=integration,
    control.compute=list(config=TRUE),num.threads='1:1',verbose=FALSE,
    working.directory=NULL,keep=FALSE)
  elapsed <- proc.time()[3]-start
  if (length(args)>2 && args[3]=='diagnostic') saveRDS(fit,file.path(root,'tmp',paste0('inla_diagnostic_',z$case_id,'.rds')))
  s <- fit$summary.lincomb.derived
  stopifnot(nrow(s)==16,all(is.finite(s$mean)),all(is.finite(s$sd)))
  posterior <- if(z$mode=='fixed')list(mean=s$mean,sd=s$sd,lower=s[['0.025quant']],upper=s[['0.975quant']],configurations=1) else integrate_inla_configurations(fit)
  elapsed<-proc.time()[3]-start
  rows[[i]] <- list(case_id=z$case_id,mode=z$mode,mean=unname(posterior$mean),sd=unname(posterior$sd),lower=unname(posterior$lower),upper=unname(posterior$upper),configurations=posterior$configurations,interpolated_mean_max_difference=max(abs(posterior$mean-s$mean)),seconds=unname(elapsed),status='ok',inla_version=as.character(packageVersion('INLA')))
  # Atomic checkpoint after every case. Restart uses a new shard, never mixes records.
  write_json(rows,output,auto_unbox=TRUE,digits=16,pretty=FALSE)
  cat(sprintf('INLA %s %d/%d %s\n',z$mode,i,length(inp),z$case_id));flush.console()
}
