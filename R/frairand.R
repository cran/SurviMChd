#' Frailty with random effects in high dimensional data with MCMC
#'
#' @description Random effects frailty model
#'
#' @param m Starting column number form where study variables to be selected.
#' @param n Ending column number till where study variables will get selected.
#' @param Inst Variable name of Institute information.
#' @param Delta Variable name containing the event information.
#' @param Time Variable name containing the time information.
#' @param T.min Variable name containing the time of event information.
#' @param chains Numner of MCMC chains.
#' @param iter Define number of iterations as number.
#' @param data High dimensional data having survival duration, event information and column of time for death cases.
#' @details By given m and n, a total of 3 variables can be selected.
#' @return frairandout
#' omeg[i] are frailty effects.
#'
#'
#' deviance is the model diagnostic criteria.
#' @import R2OpenBUGS
#' @export
#' @author Atanu Bhattacharjee and Akash Pawar
#' @examples
#' \dontrun{
#' ##
#' data(frailty)
#' frairand(5,7,Inst="institute",Delta="del",Time="timevar",T.min="time.min",2,6,frailty)
#' ##
#' }
#' @seealso fraidm fraidpm
#' @references Tawiah, R., Yau, K. K., McLachlan, G. J., Chambers, S. K., &
#' Ng, S. K. (2019). Multilevel model with random effects for clustered survival
#' data with multiple failure outcomes. \emph{Statistics in medicine}, \bold{38(6)}, 1036-1055.
#'
frairand <- function(m,n,Inst,Delta,Time,T.min,chains,iter,data){
  m1<-m
  m2=m1+1
  m3=n
  re <- data
  burn=(iter/2)

  ##Model in text format to be used in coxph univariate analysis
  frailtymodel1<-function(){
    for(i in 1:v){
      t[i] ~ dweib(gamma, mu[i])
      L[i] <- pow(f[i],delta[i])*pow(S[i],1-delta[i])
      LL[i] <- log(L[i]); invL[i] <- 1/L[i]
      S[i]  <- exp(-mu[i]*pow(t[i],gamma)); f[i]  <- mu[i]*gamma*pow(t[i],gamma-1)*S[i]
      log(mu[i]) <- b0 + b[1]*var1[i]+b[2]*var2[i]+b[3]*var3[i]+omeg[inst[i]]
      var2[i] ~ dnorm(mu.karno1,tau.karno[1])
      var3[i] ~ dnorm(mu.karno2[i],tau.karno[2])
      mu.karno2[i] <- c[1]+c[2]*var2[i]
    }
    gamma ~ dgamma(1,0.001); 	tau.omeg ~ dgamma(1,0.001); mu.karno1 ~ dflat()
    for (j in 1:2){
      tau.karno[j] ~ dgamma(1,0.001)
    }
    b0 ~ dnorm(0,0.001)
    for (j in 1:3){
      b[j] ~ dnorm(0,0.001)
    }
    for (j in 1:2){
      c[j] ~ dnorm(0,0.001)
    }
    for (j in 1:18){
      omeg[j] ~ dnorm(0,tau.omeg)
    }
  }

  ## some temporary file name:
  modelfile1 <- file.path(tempdir(), "frailtymodel1.txt")
  ## write model file:
  write.model(frailtymodel1, modelfile1)


  mdata <- list(inst=re[,Inst], delta=re[,Delta], t=re[,Time],	t.min=re[,T.min], var1=re[,m1],
                var2=re[,m2], var3=re[,m3], v=nrow(re))

  inits<-function(){list(mu.karno1=50,c=c(0,0),gamma=1,b0=-5,b=c(0,0,0),tau.omeg=1)}
  params1 <- c("c","omeg")
  v1.sim <- bugs(mdata, inits, model.file = modelfile1,
                 parameters.to.save = params1,
                 n.chains = chains, n.iter = iter)
  frairandout <- data.frame(v1.sim$summary)
  colnames(frairandout) <- c("Posterior Means","SD","2.5%","25%","50%","75%","97.5%","Rhat","n.eff")
  return(frairandout)
}
utils::globalVariables(c("b0","omeg"))

