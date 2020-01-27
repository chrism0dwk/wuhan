#' Import contact matrix
#'
#' @param filename the filename containing the matrix, csv with first column as row.names
#' @param replace_NA_w_zero replace NAs with zero
#'
#' @return a matrix
#' @export
import_matrix = function(filename, replace_NA_w_zero=FALSE, ...) {
  tmp = read.csv(filename, stringsAsFactors=FALSE, ...)
  n = ncol(tmp)-1
  M = data.matrix(tmp[,1+(1:n)])
  rownames(M) = tmp[,1]
  if (isTRUE(replace_NA_w_zero)) {
    M[is.na(M)] = 0
  }
  M
}


# Dataset documentation

#' Population sizes of cities in China
#'
#' @format A vector of length 187
#' \describe{
#'   each element corresponds to the population of a city in China
#'   mentioned in \code{china_cases}
#' }
#'
"china_population"


#' Case detections per day in China
#'
#' A dataset of detected cases of 2019-nCoV in Chinese
#' cities between 1st Jan 2020 and 24th Jan 2020.
#'
#' @format A data frame with 187 rows as 24 columns
#' \describe{
#'   columns represent dates, rows represent cities
#' }
#'
"china_cases"


#' Case detections per day in non-China countries
#'
#' A dataset of detected cases of 2019-nCoV in non-China
#' countries between 1st Jan 2020 and 24th Jan 2020.
#'
#' @format A data frame with 201 rows and 24 columns
#' \describe{
#'   columns represent dates, rows represent countries.
#' }
#'
"world_cases"


#' Bootstrap samples of estimated parameters from
#' the Lancaster Wuhan 2019-nCoV model.
#'
#' A dataset containing bootstrapped samples of the model parameters
#' \eqn{\beta}, \eqn{\gamma}, \eqn{I_{Wuhan,0}}, \eqn{phi} with the
#' latent period \eqn{\alpha=1/4}.
#'
#' @format A data frame with 1000 rows (samples) and 4 columns
#' \describe{
#'   \item{beta} samples of \eqn{log(\beta)}
#'   \item{gamma} samples of \eqn{log(\gamma)}
#'   \item{IOW} samples of \eqn{log(I_{Wuhan,0})}
#'   \item{phi} samples of \eqn{logit(phi)}
#' }
#'
"p_hat_alpha4.0"


#' Bootstrap samples of estimated parameters from
#' the Lancaster Wuhan 2019-nCoV model.
#'
#' A dataset containing bootstrapped samples of the model parameters
#' \eqn{\beta}, \eqn{\gamma}, \eqn{I_{Wuhan,0}}, \eqn{phi} with the
#' latent period \eqn{\alpha=1/4.4}.
#'
#' @format A data frame with 1000 rows (samples) and 4 columns
#' \describe{
#'   \item{beta} samples of \eqn{log(\beta)}
#'   \item{gamma} samples of \eqn{log(\gamma)}
#'   \item{IOW} samples of \eqn{log(I_{Wuhan,0})}
#'   \item{phi} samples of \eqn{logit(phi)}
#' }
#'
"p_hat_alpha4.4"


#' Bootstrap samples of estimated parameters from
#' the Lancaster Wuhan 2019-nCoV model.
#'
#' A dataset containing bootstrapped samples of the model parameters
#' \eqn{\beta}, \eqn{\gamma}, \eqn{I_{Wuhan,0}}, \eqn{phi} with the
#' latent period \eqn{\alpha=1/5}.
#'
#' @format A data frame with 1000 rows (samples) and 4 columns
#' \describe{
#'   \item{beta} samples of \eqn{log(\beta)}
#'   \item{gamma} samples of \eqn{log(\gamma)}
#'   \item{IOW} samples of \eqn{log(I_{Wuhan,0})}
#'   \item{phi} samples of \eqn{logit(phi)}
#' }
#'
"p_hat_alpha5.0"


#' Bootstrap samples of estimated parameters from
#' the Lancaster Wuhan 2019-nCoV model.
#'
#' A dataset containing bootstrapped samples of the model parameters
#' \eqn{\beta}, \eqn{\gamma}, \eqn{I_{Wuhan,0}}, \eqn{phi} with the
#' latent period \eqn{\alpha=1/6}.
#'
#' @format A data frame with 1000 rows (samples) and 4 columns
#' \describe{
#'   \item{beta} samples of \eqn{log(\beta)}
#'   \item{gamma} samples of \eqn{log(\gamma)}
#'   \item{IOW} samples of \eqn{log(I_{Wuhan,0})}
#'   \item{phi} samples of \eqn{logit(phi)}
#' }
#'
"p_hat_alpha6.0"
