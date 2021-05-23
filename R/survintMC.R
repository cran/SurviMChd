#' High dimensional survival analysis with interval censored data by MCMC
#'
#' @description Performs survival analysis with MCMC on a data set by computing survival interval given left and right censoring time.
#'
#' @details The survival columns of the data should be arranged as follows -
#' leftcensoring The column containing the left censoring information, must be named as 'Leftcensor'.
#' Rightcensor The column containing the right censoring information, must be named as 'Rightcensor' i.e. OS.
#' Death The column containing the death and alive information, must be names as 'Status'.
#
#' @param m Starting column number from where variables of high dimensional data will get selected.
#' @param n Ending column number till where variables of high dimensional data will get selected.
#' @param Leftcensor "Variable/column name" containing the left censoring information.
#' @param OS "Variable/column name" containing survival duration event observations.
#' @param Death "Variable/column name" containing the survival event information. i.e. Death
#' @param iter Number of MCMC iterations.
#' @param data High dimensional data containing the Left censoring, Right censoring, Status and DEG observations.
#'
#' @return survintMCout A table containing HR and CI for respective covariates.
#' @import ICBayes
#' @import icenReg
#' @importFrom graphics lines
#'
#' @author Atanu Bhattacharjee and Akash Pawar
#' @examples
#' \dontrun{
#' ##
#' data(hnscc)
#' survintMC(m=7,n=11,Leftcensor="leftcensoring",OS="os",Death="death",iter=6,data=hnscc)
#' ##
#' }
#' @seealso survMC
#' @references Bogaerts, K., Komarek, A., & Lesaffre, E. (2017). Survival analysis
#'  with interval-censored data: A practical approach with examples in R, SAS,
#'  and BUGS. CRC Press.
#' @export
survintMC <- function(m,n,Leftcensor=NULL,OS,Death,iter,data){

  data1 <- data
  data2 <- data
  colnames(data2) <- NULL
  hrt <- matrix(ncol = 1)
  hrci <- matrix(ncol = 2)
  variables <- matrix(ncol = 1)
  burn=(iter/2)
  lf<-data1[,Leftcensor]
  if(length(Leftcensor)==0)
    {
    lftcen<-c(rep(0,nrow(data)))
    }
  if(length(Leftcensor)==1)
    {
    lftcen<-lf
    }
  for(i in m:n){
    breastICB <- ICBayes(model = "case2ph",
                         L = lftcen, R = data1[,OS], status = data1[,Death],
                         xcov = data2[, i], x_user = c(0, 1),
                         knots = seq(0.1, 60.1, length = 4),
                         grids = seq(0.1, 60.1, by = 1),
                         niter = iter, burnin = (iter/2)
    )
    ngrid <- length(breastICB$S0_m)
    #plot(breastICB$grids, breastICB$S_m[1:ngrid], type = "l",
        # lty = "solid", xlab = "Survival times (months)", main = names(data1)[i],
        # ylab = "Estimated survival distributions", ylim = c(0, 1))
    #lines(breastICB$grids, breastICB$S_m[(ngrid+1):(2*ngrid)],
       #   lty = "dashed")
    HR <- exp(breastICB$coef)
    HR.CI <- exp(breastICB$coef_ci)

    hrt<-rbind(hrt,HR)
    hrci<-rbind(hrci,HR.CI)
    variables<-rbind(variables,names(data1)[i])
  }
  survintMCout <- data.frame(variables,hrt,hrci)
  survintMCout  <- survintMCout[-1,]
  colnames(survintMCout ) <- c("Variables","HR","LCL","UCL")
  rownames(survintMCout ) <- NULL
  return(survintMCout)
}
utils::globalVariables(c("lines"))
