#' Frailty in high dimensional survival data.
#'
#' @description Data set listing institutional wise survival outcomes
#'
#' @description Survival observations data for frailty model functions of SurviMChd
#' @usage data(frailty)
#' @format A \code{tibble} with 7 columns and 272 rows which are :
#' \describe{
#' \item{institute}{Institute of the sample observations}
#' \item{del}{Numberic values 0 or 1 containing death/event information}
#' \item{timevar}{Survival duration}
#' \item{time.min}{Minimum survival}
#' \item{female}{Covariate_1, gender variable indicating either a female or not}
#' \item{ph.karno}{Covariate_2}
#' \item{pat.karno}{Covariate_3}}
#'
#' @examples data(frailty)
#'
"frailty"
