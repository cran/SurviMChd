#' Survival analysis using Cox Proportional Hazards with MCMC.
#'
#' @description Performs survival analysis using Cox Proportional Hazards with MCMC.
#' @details The survival columns of the data should be arranged as follows -
#' Death Death status=1 if died otherwise 0.
#' OS Survival duration measured as 'OS'
#' t.len Number of censored times
#'
#' @param m Starting column number from where variables of high dimensional data will get selected.
#' @param n Ending column number till where variables of high dimensional data will get selected.
#' @param Time Variable/Column name containing the information on duration of survival
#' @param Event Variable/Column name containing the information of survival event
#' @param chains Number of chains to perform
#' @param adapt Number of adaptations to perform
#' @param iter Number of iterations to perform
#' @param data High dimensional data having survival duration and event.
#' @return Data set containing Posterior HR estimates, SD and quantiles.
#' @import rjags
#'
#' @author Atanu Bhattacharjee and Akash Pawar
#' @examples
#' \donttest{
#' ##
#' data(mcsurv)
#' survMC(m=4,n=8,Time="OS",Event="Death",chains=2,adapt=100,iter=1000,data=mcsurv)
#' ##
#' }
#' @seealso survintMC
#' @references Bhattacharjee, A. (2020). Bayesian Approaches in Oncology Using R
#' and OpenBUGS. CRC Press.
#' @export
survMC <- function(m,n,Time,Event,chains,adapt,iter,data)
{
  if(Time!="OS"){
    names(data)[names(data) == Time] <- "OS"
  }
  if(Event!="Death"){
    names(data)[names(data) == Event] <- "Death"
  }


  data<-data[order(data$OS),]
  var1 <- colnames(data)
  nr <- nrow(data)

  data1 <- subset(data, Death == 1) #subsetting data with death status = 1
  u <- unique(data1$OS) #creating a vector with unique values of OS

  #adding a condition for censoring time vector to include the last censored patient when censoring = 0

  if ((data$Death[nrow(data)])==0){
    u1<-c(u,data$OS[nrow(data)])
  } else {
    u1 <- u
  }
  u2 <- sort(u1)
  t.len<-(length(u2)-1)



  model_jags <- "
  data{
    # Set up data
  for(i in 1:N) {
    for(j in 1:T) {
    Y[i,j] <- step(obs.t[i] - t[j] + eps)
    dN[i, j] <- Y[i, j] * step(t[j + 1] - obs.t[i] - eps) * fail[i]
    }
  }
  }

  # Model
  model{
  for(i in 1:N){
    betax[i,1] <- 0
    for(k in 2:(p+1)){
      betax[i,k] <- betax[i,k-1] + beta[k-1]*x[i,k-1]
    }
  }
  for(j in 1:T) {
    for(i in 1:N) {
    dN[i, j] ~ dpois(Idt[i, j]) # Likelihood
    Idt[i, j] <- Y[i, j] * exp(betax[i,p+1]) * dL0[j] # Intensity
    }
    dL0[j] ~ dgamma(mu[j], c)
    mu[j] <- dL0.star[j] * c # prior mean hazard
  }
  c <- 0.001
  r <- 0.1
  for (j in 1 : T) {
    dL0.star[j] <- r * (t[j + 1] - t[j])
  }
  for(k in 1:p){
    beta[k] ~ dnorm(0.0,0.000001)
  }
  }"

  params <- c("beta","dL0")

  inits <-  function(){list( beta = rep(0,p), dL0 = rep(0.0001,bigt))}

  x2 <- rep(0,nrow(data))

  q <- matrix(nrow=0,ncol=5)
  s <- matrix(nrow=0,ncol=2)
  di <- matrix(nrow=0,ncol=1)
  for(i in m:n){
    x1 <- data[(1:nrow(data)),i]
    x = t(rbind(x1,x2))

    datafi <- list(x=x,obs.t=data$OS,t=u2,T=t.len,N=nrow(data),fail=data$Death,eps=1E-10,p=2)

    jags <- jags.model(textConnection(model_jags),
                       data = datafi,
                       n.chains = chains,
                       n.adapt = adapt)


    samps <- coda.samples(jags, params, n.iter=iter)
    s1 <- summary(samps)
    stats <- s1$statistics[1,c(1:2)]
    s <- rbind(s,stats)
    quan <- s1$quantiles[1,]
    q <- rbind(q,quan)
    d = dic.samples(jags, n.iter=iter)
    meandeviance <- round(sum(d$deviance),2)
    di <- rbind(di,meandeviance)
  }
  results <- cbind(s,q)
  expresults <- exp(results)

  Variables <- names(data)[m:n]
  expresults <- data.frame(Variables,expresults,di)
  colnames(expresults)<-c("Variable","Posterior Means","SD","2.5%","25%","50%","75%","97.5%","DIC")
  rownames(expresults) <- NULL
  return(expresults)
}
utils::globalVariables(c("Death","N","step","obs.t","eps","fail","x1","dL0","p","bigt"))
