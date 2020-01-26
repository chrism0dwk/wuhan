#' Returns the logit of x
#'
#' @param x input value
#'
#' @return logit of x
logit = function(x) {
  log(x/(1-x))
}

#' Returns the inverse logit of x
#'
#' @param x input value
#'
#' @return inverse logit of x
invlogit = function(x) {
  expx = exp(x)
  expx / (1 + exp(x))
}


#' Likelihood function for Wuhan nCoV-2019
#'
#' This function returns a function encoding the likelihood for the nCoV-2019
#' model.  It assumes Poisson-distributed increments of case
#' reports in China, and Poisson-distributed increments in numbers of
#' infected passengers on planes elsewhere.
#'
#' @param y a \eqn{n \times T} matrix of case reports in China
#' @param z a \eqn{m \times T} matrix of case reports elsewhere
#' @param N the population sizes within China, length n
#' @param K the within-China air travel matrix, \eqn{n \times n}
#' @param W the international air travel matrix, \eqn{n \times m}
#' @param sim_fun a function which returns a simulation from a disease model
#' @param phi_mask a vector of 0s and 1s determining which cities in China to apply the
#'     underreporting parameter \eqn{phi} to.
#' @param agg_up_to aggregate the first however many case detection records, allowing for
#'     a delay in receiving counts from the first few cases.
#'
#' @details This function returns a closure -- another function that encapsulates the
#' data passed to the containing function.  See the return value for the function signature.
#'
#' @return a function to calculate the log likelihood.  This function has signature
#' \code{logp_fn(params, visualise=FALSE)} where \code{params} is a vector of parameters
#' \code{c(beta, gamma, I0W, phi)}.  If \code{visualise} is \code{TRUE}, then parameter values
#' are printed to the console and a graph showing how the ODE mean function matches the observed
#' timeseries (in Wuhan) is displayed.  This is useful for tracking the progression of various
#' optimisers.
#'
#' @export
#' @import assertthat
LogLikelihood = function(y, z, N, K, W, sim_fun,
                         phi_mask=(rownames(K)=='Wuhan'), agg_up_to=11) {

  # Transpose both y and z for consistency with ODE output
  y = t(y)
  z = t(z)

  # Calculates the log likelihood
  #
  # @param param a vector of model parameters c(beta, gamma, I0W, phi)
  # @param visualise if TRUE then print out model parameters and a graph of
  # modelled case detections in Wuhan.  Useful for following progress of optimisers
  #
  # @return the log likelihood of the data conditional on the model and parameters
  #
  logp_fn = function(param, visualise=FALSE) {
    sim_param = exp(param[1:3])

    if(isTRUE(visualise)) {
      pparam = c(beta=exp(param[1]), gamma=exp(param[2]),
                 I0W=exp(param[3]), phi=invlogit(param[4]))
      print(pparam)
    }

    expected = sim_fun(sim_param)
    p_detect = rep(1, length(N))
    p_detect[phi_mask] = invlogit(param[4]) # phi, underreporting

    # Observe increments in R in China
    exp_incr = t(t(diff(expected$R)) * p_detect)

    # China observation model
    y_prime = y[agg_up_to:nrow(y),]
    exp_incr_prime = rbind(colSums(exp_incr[1:agg_up_to,]),
                           exp_incr[(agg_up_to+1):nrow(exp_incr),])
    assertthat::assert_that(all(dim(y_prime)==dim(exp_incr_prime)))
    llik_china = dpois(y_prime, exp_incr_prime, log=TRUE)
    llik_china = sum(llik_china)

    if (isTRUE(visualise)) {
      plot(cumsum(y_prime[, phi_mask]))
      lines(cumsum(exp_incr_prime[, phi_mask]), col=2)
    }

    # Rest of world observation model
    china_prev = t(t(expected$I / N) * p_detect)
    flight_prev = china_prev %*% W
    llik_world = sum(dpois(z, flight_prev, log=TRUE))

    llik_china + llik_world
  }
}

