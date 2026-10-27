# Integrate actual INLA Gaussian configurations rather than interpolated
# univariate density summaries. Weights and conditional means/precisions come
# from R-INLA; the independent Python quadrature is used only for validation.
integrate_inla_configurations <- function(fit) {
  configs <- fit$misc$configs
  joint <- as.matrix(fit$joint.hyper)
  stopifnot(ncol(joint)==3,configs$ntheta==1)
  theta <- vapply(configs$config,function(z) as.numeric(z$theta),numeric(1))
  weights <- joint[,3];weights<-weights/sum(weights)
  stopifnot(all(is.finite(weights)),all(weights>=0))
  cont <- configs$contents
  offset <- sum(cont$length)-configs$n
  first <- setNames(cont$start-offset,cont$tag)
  mu <- first['mu']+0:7;delta<-first['delta']+0:15
  C <- matrix(0,16,configs$n)
  for (j in 1:8) {C[j,mu[j]]<-1;C[j,delta[j]]<-1;C[j+8,mu[j]]<-1;C[j+8,delta[j+8]]<-1}
  means<-variances<-matrix(NA_real_,nrow(joint),16)
  for (h in seq_len(nrow(joint))) {
    pick<-which.min(abs(theta-joint[h,1]));stopifnot(abs(theta[pick]-joint[h,1])<1e-6)
    z<-configs$config[[pick]]
    q<-as.matrix(z$Q)
    # INLA stores one triangle of the latent precision.
    upper<-sum(abs(q[upper.tri(q)]));lower<-sum(abs(q[lower.tri(q)]))
    if (!isSymmetric(q,tol=1e-9)) {
      stopifnot(min(upper,lower)<1e-8)
      q<-as.matrix(forceSymmetric(z$Q,uplo=if(upper>lower)'U' else 'L'))
    }
    means[h,]<-as.vector(C%*%z$mean)
    variances[h,]<-diag(C%*%solve(q,t(C)))
  }
  stopifnot(all(is.finite(means)),all(variances>0))
  mean<-as.vector(crossprod(weights,means))
  sd<-sqrt(pmax(0,as.vector(crossprod(weights,variances+means^2))-mean^2))
  limits<-matrix(NA_real_,16,2)
  for(j in 1:16) {
    bound<-max(1,abs(means[,j])+12*sqrt(variances[,j]))
    for(k in 1:2) limits[j,k]<-uniroot(function(x)sum(weights*pnorm((x-means[,j])/sqrt(variances[,j])))-c(.025,.975)[k],c(-bound,bound),tol=1e-10)$root
  }
  list(mean=mean,sd=sd,lower=limits[,1],upper=limits[,2],configurations=nrow(joint))
}
