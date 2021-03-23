#' Frailty with Discrete Mixture  Model
#'
#' @description Discrete mixture model with MCMC
#'
#' @param m Starting column number form where study variables to be selected.
#' @param n Ending column number till where study variables will get selected.
#' @param Inst Variable name of Institute information.
#' @param Delta Variable name containing the event information.
#' @param Time Variable name containing the time information.
#' @param T.min Variable name containing the time of event information.
#' @param chains Number of MCMC chains
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
#' deviance is the model diagnostic criteria.
#'
#' @import R2OpenBUGS
#' @export
#' @examples
#' \dontrun{
#' ##
#' data(frailty)
#' fraidm(5,7,Inst="institute",Delta="del",Time="timevar",T.min="time.min",2,7,frailty)
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
fraidm <- function(m,n,Inst,Delta,Time,T.min,chains,iter,data){
  m1<-m
  m2=m1+1
  m3=n
  re=data
  burn=(iter/2)

  frailtymod2 <- function(){

    for (i in 1 : v) {t[i] ~ dweib(gamma, mu[i])
      L[i] <- pow(f[i],delta[i])*pow(S[i],1-delta[i]);   LL[i] <- log(L[i]); invL[i] <- 1/L[i]
      S[i]  <- exp(-mu[i]*pow(t[i],gamma));      f[i]  <- mu[i]*gamma*pow(t[i],gamma-1)*S[i]
      var2[i] ~ dnorm(mu.karno1,tau.karno[1])
      var3[i] ~ dnorm(mu.karno2[i],tau.karno[2])
      mu.karno2[i] <- c[1]+c[2]*var2[i]
      log(mu[i]) <- b[1]*var1[i]+b[2]*var2[i]+b[3]*var3[i]+omeg[G[inst[i]]]}
    gamma ~ dgamma(1,0.001);
    for (j in 1:2) {tau.karno[j] ~ dgamma(1,0.001)}
    omeg[1] ~ dnorm(0,0.001); omeg[2] <- omeg[1]+del
    del ~ dexp(1)
    mu.karno1 ~ dflat()
    for (j in 1:3) {b[j] ~ dnorm(0,0.001)}
    for (j in 1:2) {c[j] ~ dnorm(0,0.001)}
    for (j in 1:18) {G[j] ~ dcat(pi[1:2]); Gcat[j] <- equals(G[j],1)}
    pi[1:2] ~ ddirch(alph[1:2])
    for (j in 1:2) {alph[j] <- 0.5}

  }

  ## some temporary file name:
  modelfile2 <- file.path(tempdir(), "frailtymod2.txt")
  ## write model file:
  write.model(frailtymod2, modelfile2)


  mdata <- list(
    inst=re[,Inst], delta=re[,Delta], t=re[,Time],	t.min=re[,T.min],var1=re[,m1],
    var2=re[,m2], var3=re[,m3], v=nrow(re))

  inits<-function(){list(mu.karno1=50,c=c(0,0),gamma=1,b=c(0,0,0),omeg=c(0,NA))}
  params1 <- c("b","c","omeg")
  v1.sim <- bugs(mdata, inits, model.file = modelfile2,
                 parameters.to.save = params1,
                 n.chains = chains, n.iter = iter)

  fraidmout <- data.frame(v1.sim$summary)
  colnames(fraidmout) <- c("Posterior Means","SD","2.5%","25%","50%","75%","97.5%","Rhat","n.eff")
  return(fraidmout)
}
utils::globalVariables(c("v","pow","delta","var2","log<-","b","var1","var3","G","inst","del","equals"))
