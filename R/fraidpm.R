#' Frailty with drichlet process mixture
#'
#' @description Frailty analysis on high dimensional data by Drichlet process mixture.
#'
#' @param m Starting column number form where study variables to be selected.
#' @param n Ending column number till where study variables will get selected.
#' @param Inst Variable name of Institute information.
#' @param Delta Variable name containing the event information.
#' @param Time Variable name containing the time information.
#' @param T.min Variable name containing the time of event information.
#' @param chains Number of MCMC chains.
#' @param iter Define number of iterations as number.
#' @param data High dimensional data, event information given as (delta=0 if alive, delta=1 if died). If patient is censored then t.min=duration of survival. If patient is died then t.min=0. If patient is died then t=duration of survival. If patient is alive then t=NA.
#'
#' @details By given m and n, a total of 3 variables can be selected.
#' @return fraidpmout
#' omeg[i] are frailty effects.
#'
#'
#' deviance is the model diagnostic criteria.
#'
#' @import R2OpenBUGS
#' @importFrom stats cor step
#' @export
#' @author Atanu Bhattacharjee and Akash Pawar
#' @examples
#' \dontrun{
#' ##
#' data(frailty)
#' fraidpm(5,7,Inst="institute",Delta="del",Time="timevar",T.min="time.min",2,6,frailty)
#' ##
#' }
#' @seealso fraidm frairand
#' @references Bhattacharjee, A. (2020). Bayesian Approaches in Oncology Using
#' R and OpenBUGS. CRC Press.
#' @references Congdon, P. (2014). Applied bayesian modelling (Vol. 595).
#' John Wiley & Sons.
#'
fraidpm <- function(m,n,Inst,Delta,Time,T.min,chains,iter,data){
  m1<-m
  m2=m1+1
  m3=n
  re <- data
  burn = (iter/2)
  #openbugs model function
  frailtymod3 <- function(){
    for (i in 1 : N) {t[i] ~ dweib(gamma, mu[i])
      L[i] <- pow(f[i],delta[i])*pow(S[i],1-delta[i]);   LL[i] <- log(L[i]); invL[i] <- 1/L[i]
      S[i]  <- exp(-mu[i]*pow(t[i],gamma));      f[i]  <- mu[i]*gamma*pow(t[i],gamma-1)*S[i]
      var2[i] ~ dnorm(mu.karno1,tau.karno[1])
      var3[i] ~ dnorm(mu.karno2[i],tau.karno[2])
      mu.karno2[i] <- c[1]+c[2]*var2[i]
      log(mu[i]) <- b[1]*var1[i]+b[2]*var2[i]+b[3]*var3[i]+omeg[inst[i]]}
    gamma ~ dgamma(1,0.001);
    for (j in 1:2) {tau.karno[j] ~ dgamma(1,0.001)}
    mu.karno1 ~ dflat();     for (j in 1:2) {c[j] ~ dnorm(0,0.001)}
    for (j in 1:3) {b[j] ~ dnorm(0,0.001)}
    # Cluster Allocation under DPP
    for (j in 1:K) {G[j] ~ dcat(p[1:M])
      # frailty in institution j
      omeg[j] ~ dnorm(nu[G[j]],tau.omeg[G[j]])
      for (k in 1:M) {post[j,k] <- equals(G[j],k)}}
    # baseline prior
    nu0 ~ dnorm(0,0.001);  tau0 ~ dgamma(1,0.001)
    for (i in 1:M) {nu[i] ~ dnorm(nu0,tau0)
      tau.omeg[i] ~ dgamma(1,0.001); }
    # truncated Dirichlet process
    alpha <- 5;          V[M] <- 1
    for (k in 1:M-1){  V[k] ~ dbeta(1,alpha)}
    p[1] <- V[1]
    for (j in 2:M)    {  p[j] <- V[j]*(1-V[j-1])*p[j-1]/V[j-1]}
    # total clusters
    TClus <- sum(NonEmp[])
    for (j in 1:M) {NonEmp[j] <- step(sum(post[,j])-1)}
  }


  ## some temporary filename:
  modelfile3 <- file.path(tempdir(), "frailtymod3.txt")
  ## write model file:
  write.model(frailtymod3, modelfile3)

  mdata <- list(N=nrow(re),M=25,K=18,
                inst=re[,Inst], delta=re[,Delta], t=re[,Time],	t.min=re[,T.min], var1=re[,m1],
                var2=re[,m2], var3=re[,m3])

  inits<-function(){list(mu.karno1=50,c=c(0,0),gamma=1,b=c(0,0,0),nu0=-5,tau0=1)}

  params1 <- c("omeg")
  v1.sim <- bugs(mdata, inits, model.file = modelfile3,
                 parameters.to.save = params1,
                 n.chains = chains, n.iter = iter)

  fraidpmout <- data.frame(v1.sim$summary)
  colnames(fraidpmout) <- c("Posterior Means","SD","2.5%","25%","50%","75%","97.5%","Rhat","n.eff")
  return(fraidpmout)
}
utils::globalVariables(c("N","omeg","K","M","step"))

