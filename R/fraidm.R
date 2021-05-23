#' Frailty with Discrete Mixture  Model
#'
#' @description Discrete mixture model with MCMC
#'
#' @param m Starting column number form where study variables to be selected.
#' @param n Ending column number till where study variables will get selected.
#' @param Ins Variable name of Institute information.
#' @param Del  Variable name containing the event information.
#' @param Time Variable name containing the time information.
#' @param T.min Variable name containing the time of event information.
#' @param chn Number of MCMC chains
#' @param iter Define number of iterations as number.
#' @param data High dimensional data, event information given as (delta=0 if alive, delta=1 if died). If patient is censored then t.min=duration of survival. If patient is died then t.min=0. If patient is died then t=duration of survival. If patient is alive then t=NA.
#'
#'
#' @details By given m and n, a total of 3 variables can be selected.
#'
#' @return fraidmout - b[1] is the posterior estimate of the regression coefficient for first covariate.
#'
#'
#' b[2] is the posterior estimate of the regression coefficient for second covariate.
#'
#'
#' b[3] is the posterior estimate of the regression coefficient for third covariate.
#'
#'
#' omega[1] and omega[2] are frailty effects.
#'
#'
#' c[1] and c[2] are regression intercept and coefficients of covariates over mean effect.
#'
#'
#' @import rjags
#' @export
#' @examples
#' \dontrun{
#' ##
#' data(frailty)
#' fraidm(m=5,n=7,Ins="institute",Del="del",Time="timevar",T.min="time.min",chn=2,iter=6,data=frailty)
#' ##
#' }
#' @seealso fraidpm frairand
#' @references Bhattacharjee, A. (2020). Bayesian Approaches in Oncology Using
#' R and OpenBUGS. CRC Press.
#' @references Congdon, P. (2014). Applied bayesian modelling (Vol. 595).
#' John Wiley & Sons.
#'
#' @references {
#' }
fraidm <- function(m,n,Ins,Del,Time,T.min,chn,iter,data){
  Inst<-Ins
  Delta<-Del
  chains<-chn
  m1<-m
  m2=m1+1
  m3=n
  re=data
  mdata <- list(
    inst=re[,Inst], delta=re[,Delta], t=re[,Time],	#t.min=re[,T.min],
    var1=re[,m1],
    var2=re[,m2], var3=re[,m3], v=nrow(re))

  dat1 <- mdata


  f1 <- "model{
    for (i in 82 : v) {t[i] ~ dweib(gamma, mu[i])
      L[i] <- pow(f[i],delta[i])*pow(S[i],1-delta[i])
      LL[i] <- log(L[i]); invL[i] <- 1/L[i]
      S[i]  <- exp(-mu[i]*pow(t[i],gamma))
      f[i]  <- mu[i]*gamma*pow(t[i],gamma-1)*S[i]
      var2[i] ~ dnorm(mu.karno1,tau.karno[1])
      var3[i] ~ dnorm(mu.karno2[i],tau.karno[2])
      mu.karno2[i] <- c[1]+c[2]*var2[i]
      log(mu[i]) <- b[1]*var1[i]+b[2]*var2[i]+b[3]*var3[i]+omeg[G[inst[i]]]}
    gamma ~ dgamma(1,0.001);
    for (j in 1:2) {tau.karno[j] ~ dgamma(1,0.001)}
    omeg[1] ~ dnorm(0,0.001); omeg[2] <- omeg[1]+del
    del ~ dexp(1)
    mu.karno1 ~ dnorm(0,1)
    for (j in 1:3) {b[j] ~ dnorm(0,0.001)}
    for (j in 1:2) {c[j] ~ dnorm(0,0.001)}
    for (j in 1:18) {G[j] ~ dcat(pi[1:2]); Gcat[j] <- equals(G[j],1)}
    pi[1:2] ~ ddirch(alph[1:2])
    for (j in 1:2) {alph[j] <- 0.5}
  }"

  inits1<-list(mu.karno1=50,c=c(0,0),gamma=1,b=c(0,0,0),omeg=c(0,NA))

  dat1 <- mdata
  jagsmod2 <- jags.model(textConnection(f1),
                         data = dat1,
                         inits = inits1,
                         n.chains = chains)
  params1 <- c("b","c","omeg")
  samps <- coda.samples(jagsmod2, params1, n.iter=iter)
  s1 <- summary(samps)
  stats <- data.frame(s1[1]$statistics[,c(1,2)])
  quant <- data.frame(s1[2]$quantiles[,])
  statsquant <- data.frame(stats,quant)
  return(statsquant)
}
