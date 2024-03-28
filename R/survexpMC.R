#' Exponential survival analysis with MCMC
#'
#' @description Survival analysis with exponential distribution by MCMC
#'
#' @param m1 Starting column number from where variables of high dimensional data will be selected.
#' @param n1 Ending column number till where variables of high dimensional data will get selected.
#' @param m2 Starting column number from where demographic observations starts
#' @param n2 Ending column number of the demographic observations
#' @param chains Number of MCMC chains
#' @param iter Number of MCMC iterations
#' @param data High dimensional data having survival duration as (OS), event information as Death (1 if died, or 0 if alive).
#' @return survexpMCout A data set listing estimated posterior means and deviances
#' @import R2jags
#' @import dplyr
#' @importFrom stats cor step
#'
#' @author Atanu Bhattacharjee and Akash Pawar
#' @examples
#' \donttest{
#' ##
#' data(headnneck)
#' survexpMC(m1=8,n1=12,m2=4,n2=7,chains=2,iter=10,data=headnneck)
#' ##
#' }
#' @seealso survweibMC
#' @references Kumar, M., Sonker, P. K., Saroj, A., Jain, A., Bhattacharjee, A.,
#'  & Saroj, R. K. (2020). Parametric survival analysis using R: Illustration
#'  with lung cancer data. \emph{Cancer Reports}, \bold{3(4)}, e1210.
#' @export
survexpMC <- function(m1,n1,m2,n2,chains,iter,data){

  p=m1;q=n1;r=m2;s=n2;
  vr <- c(r:s)
  y1 <- subset(data,select=c(p:q))
  rownames(y1)<-NULL
  colnames(y1)<-NULL
  y1<-as.matrix(y1)
  variables = 4
  m=length(p:q)
  data <- data %>%
    mutate(surt = if_else(data$Death == 1, data$OS, 0))
  data <- data %>%
    mutate(surt.cen = if_else(data$Death == 0, data$OS, 0))
  data$surt[data$surt == 0] <- NA
  M=length(p:q)

  surt = data$"surt"
  surt.cen=data$"surt.cen"
  colnames(data)<-NULL
  randgrp1=c(data[,vr[1]])
  gender1=c(data[,vr[2]])
  stratum1=c(data[,vr[3]])
  prevoi1=c(data[,vr[4]])
  a<-list(N=nrow(data),
          M=length(p:q),
          betamu1=c(rep(0,6)),
          betamu2=c(rep(0,M)),
          Sigma1=diag(0.01,6,6),
          Sigma2=diag(0.01,M,M),
          U0=c(rep(0,2)),
          #R=diag(0.01,2,2),
          t=c(0,2,seq(6,(6*(M-2)),6)),
          Y=y1,
          surt=surt,
          #surt.cen=surt.cen,
          randgrp1=randgrp1,
          gender1=gender1,
          stratum1=stratum1,
          prevoi1=prevoi1)

  fx <- function()
      {
     for (i in 1:N) {
        for (j in 1:M) {
          Y[i, j] ~ dnorm(muy[i, j], tauz)
          muy[i, j]<-beta1[1]+
            beta1[2]*t[j]+
            beta1[3]*t[j]*randgrp1[i]+
            beta1[4]*gender1[i]+
            beta1[5]*prevoi1[i]+
            beta1[6]*stratum1[i]+
            U[i,1]+
            U[i,2]*t[j]

          yp1[i,j] ~ dnorm(muy[i, j], tauz)
          r12[i,j]<-Y[i,j]-yp1[i,j]
          sqr1[i,j]<-r12[i,j]*r12[i,j]
        }

        surt[i] ~ dweib(p,mut[i])
        log(mut[i])<-beta2[1]+beta2[2]*randgrp1[i]+beta2[3]*gender1[i]+
          beta2[4]*prevoi1[i]+beta2[5]*stratum1[i]+r1*U[i, 1]+r2*U[i, 2]

        U[i,1:2] ~ dmnorm(U0[],tau[,])
        surP[i]~ dweib(p,mut[i])
        surR[i]<-surt[i]-surP[i]
        sqr2[i]<-surR[i]*surR[i]
      }
      mspe1<-mean(sqr1[,])
      mspe2<-mean(sqr2[])
      p <- 1                  #  Use this for Exponential model
      #p ~ dgamma(1,1)  #  Use this for full Weibull model

      sigmaz<-1/tauz

      tau[1:2,1:2]<-inverse(sigma[,])

      #priors

      sigmab1 ~ dunif(0,100)
      sigmab2 ~ dunif(0,100)
      cor ~ dunif(-1,1)

      sigma[1,1] <- pow(sigmab1,2)
      sigma[2,2] <- pow(sigmab2,2)
      sigma[1,2] <- cor*sigmab1*sigmab2
      sigma[2,1] <- sigma[1,2]

      beta1[1:6]~dmnorm(betamu1[],Sigma1[,])
      tauz~dgamma(0.1, 0.1)
      beta2[1:5]~dmnorm(betamu2[],Sigma2[,])
      r1~dnorm(0, 0.01)
      r2~dnorm(0, 0.01)
  }

  #-------------------------------------------------------------------------------
  inits<-function(){list(beta1=c(rep(0,6)),tauz=1,
              beta2=c(rep(0,M)),r1=0,r2=0,
              sigmab1=1,sigmab2=1,cor=0,
              U=matrix(rep(0,(nrow(data)*2)),nrow(data),2))}

  params1 <- c("beta1")

  #### Set up the JAGS model and settings
  jags.e <- jags(data = a, inits = inits, model.file = fx, parameters.to.save = params1,
                 n.chains=chains,n.iter = iter, n.burnin = iter/2 )


  return(jags.e)

}

utils::globalVariables(c("beta1","U","yp1","Y","log<-","beta2","r1","r2","surP","tauz","inverse","pow","sigmab1","sigmab2"))

