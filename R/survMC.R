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
#' @param chains Number of chains to be considered
#' @param iter Number of iterations to be applied
#' @param data High dimensional data having survival duration and event.
#' @return survMCout Data set containing Posterior HR estimates and deviance values
#' @importFrom R2OpenBUGS bugs
#' @importFrom R2OpenBUGS write.model
#'
#' @author Atanu Bhattacharjee and Akash Pawar
#' @examples
#' \dontrun{
#' ##
#' data(mcsurv)
#' survMC(4,8,2,10,mcsurv)
#' ##
#' }
#' @seealso survintMC
#' @references Bhattacharjee, A. (2020). Bayesian Approaches in Oncology Using R
#' and OpenBUGS. CRC Press.
#' @export
survMC <- function(m,n,chains,iter,data){
  rea <- data

  re<-rea[order(rea$OS),]
  var1 <- colnames(re)
  nr <- nrow(re)
  re1 <- subset(re, Death == 1) #subsetting data with death status = 1
  u <- unique(re1$OS) #creating a vector with unique values of OS

  #adding a condition for censoring time vector to include the last censored patient when censoring = 0

  if ((re$Death[nrow(re)])==0){
    u1<-c(u,re$OS[nrow(re)])
  } else {
    u1 <- u
  }
  u2 <- sort(u1)
  t.len<-(length(u2)-1)


  #Creating dummy matrices for storing results
  summ<-matrix(nrow=0,ncol=9)
  dic<-matrix(nrow=0,ncol=1)
  vname<-matrix(nrow=0,ncol=1)

  #creating a list to store graphical results
  r=list()
  #-------------------------------------------------------------------------------
  ##Model in text format to be used in coxph univariate analysis
  uniCoxModel<-function(){
    # Set up data
    for(i in 1:N) {
      for(j in 1:T) {
        # risk set = 1 if obs.t > = t
        Y[i,j] <- step(obs.t[i] - t[j] + eps)
        # counting process jump = 1 if obs.t in [t[j], t[j+1])
        # i.e. if t[j] < = obs.t < t[j+1]
        dN[i, j] <- Y[i, j] * step(t[j + 1] -
                                     obs.t[i] - eps) * fail[i]
      }
    }
    # Model
    for(j in 1:T) {
      for(i in 1:N) {
        dN[i, j] ~ dpois(Idt[i, j]) # Likelihood
        Idt[i, j] <- Y[i, j] * exp(beta * x1[i]) *
          dL0[j] # Intensity
      }
      dL0[j]~dgamma(mu[j], c)
      mu[j] <- dL0.star[j] * c # prior mean hazard
      # Survivor function = exp(-Integral{l0(u)du})^exp(beta*z)
      S.treat[j] <- pow(exp(-sum(dL0[1 : j])),
                        exp(beta *-0.5));
      S.placebo[j] <- pow(exp(-sum(dL0[1 : j])),
                          exp(beta * 0.5));
    }
    c <- 0.001
    r <- 0.1
    for (j in 1 : T) {
      dL0.star[j] <- r * (t[j + 1] - t[j])
    }
    beta ~ dnorm(0.0,0.01)
    # hazard ratio for group 1 versus group 2
    HR<-exp(beta)
  }

  ## some temporary filename:
  modelfile <- file.path(tempdir(), "uniCoxModel.txt")
  ## write model file:
  write.model(uniCoxModel, modelfile)
  #-------------------------------------------------------------------------------
  # initial values
  inits<-function(){list(beta = 0.0,
                         dL0 = c(rep(1.0,t.len)))}
  #vector of variable names
  var2 <- var1[m:n]
  for (i in (m:n)){
    mdata1 <- list(N=nrow(re), T=t.len, eps=1.0E-10, obs.t=re$OS,
                   fail=re$Death, x1=re[(1:nrow(re)),i], t=u2)

    params1 = c("HR")
    #R2OpenBugs
    v1.sim <- bugs(mdata1, inits, model.file = modelfile,
                   parameters.to.save =  params1,
                   n.chains = chains, n.iter = iter)

    #appending summary results for each variables
    summ<-rbind(summ,v1.sim$summary)

    #appending DIC values for each variables
    dic<-rbind(dic,v1.sim$DIC)
    dic<-rbind(dic,v1.sim$DIC)

    #appending variable names
    vname<-rbind(vname,colnames(re)[i])
    vname<-rbind(vname,colnames(re)[i])


    #appending plots for each variables
    r[[colnames(re)[i]]]<-plot(v1.sim,main =" var1[i]")

  }


  #-------------------------------------------------------------------------------
  #combining results to form a data set
  dicdf<-data.frame(dic)
  vnamedf<-data.frame(vname)
  summdf<-data.frame(summ)
  survMCout <- data.frame(vname,dic,summ)

  #calling the results
  colnames(survMCout) <- c( "Variables", "DIC", "Posterior Mean", "SD",
                         "2.5%","25%","50%","75%","97.5%","Rhat","n.eff")


  return(survMCout)

}
utils::globalVariables(c("Death","N","step","obs.t","eps","fail","x1","dL0"))


