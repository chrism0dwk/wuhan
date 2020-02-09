#' Aggregate by province
#'
#' Aggregates case counts by province
#'
#' @param city_counts a \eqn{T\times M} matrix of case counts where T represents time and M represents cities.  Must have colnames denoting cities.
#' @param provinces a vector of provinces corresponding to columns in \code{city_counts}
#'
#' @return a \eqn{T\times Q} matrix of case counts where T represents time and Q represents province
#' @import reshape2
aggregate_by_province = function(city_counts, provinces)
{
  city_counts = as.data.frame(city_counts)
  cities = colnames(city_counts)
  city_counts$date = 1:nrow(city_counts)
  counts_long = reshape2::melt(city_counts, id.vars="date", variable.name="city")
  counts_long$province = provinces[match(counts_long$city, cities)]
  province_counts = group_by(counts_long, date, province) %>% summarise(cases=sum(value))
  reshape2::acast(province_counts, date~province, value.var="cases")
}


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
LogLikelihood = function(y, z, pop_info, K_early, K_late, W_early, W_late, sim_fun, t_agg=23,
                         phi_mask=(pop_info$City=='Wuhan'), first_agg=11) {

  # Transpose both y and z for consistency with ODE output
  y = t(y)
  z = t(z)
  N = pop_info$Population

  y_early = y[1:(t_agg-1),]
  y_late = aggregate_by_province(y[23:nrow(y),], pop_info$Province)

  # Calculates the log likelihood
  #
  # @param param a vector of model parameters c(beta_w_early, beta_w_late, beta_china, gamma, I0W, phi)
  # @param visualise if TRUE then print out model parameters and a graph of
  # modelled case detections in Wuhan.  Useful for following progress of optimisers
  #
  # @return the log likelihood of the data conditional on the model and parameters
  #
  logp_fn = function(param, visualise=FALSE) {
    sim_param = exp(param[1:5])
    phi = invlogit(param[6]) # Phi

    if(isTRUE(visualise)) {
      pparam = c(beta_w_early=sim_param[1], beta_w_late=sim_param[2],
                 beta_china=sim_param[3], gamma=sim_param[4],
                 I0W=sim_param[5], phi=phi)
      print(pparam)
    }

    expected = sim_fun(sim_param)
    p_detect = rep(1, length(N))
    p_detect[phi_mask] = phi # phi, underreporting

    # Observe increments in R in China
    exp_incr = t(t(diff(expected$R)) * p_detect)

    llik = numeric(4)
    # Early observations - China
    y_prime = y_early[first_agg:nrow(y_early),]
    exp_incr_prime = rbind(colSums(exp_incr[1:first_agg,]),
                           exp_incr[(first_agg+1):(t_agg-1),])
    assertthat::assert_that(all(dim(y_prime)==dim(exp_incr_prime)))
    llik_china = dpois(y_prime, exp_incr_prime, log=TRUE)
    llik[1] = sum(llik_china)

    # Early observations - Rest of world
    china_prev = t(t(expected$I / N))[-1,] # Omit starting value
    flight_prev_early = china_prev[1:(t_agg-1),] %*% W_early
    llik[2] = sum(dpois(z[1:(t_agg-1),], flight_prev_early, log=TRUE))

    # Late observations - China aggregate by province
    exp_incr_prov = aggregate_by_province(exp_incr[t_agg:nrow(exp_incr),],
                                          pop_info$Province)
    assertthat::assert_that(all(dim(y_late)==dim(exp_incr_prov)))
    llik[3] = sum(dpois(y_late, exp_incr_prov, log=TRUE))

    # Late observations - World
    flight_prev_late = china_prev[t_agg:nrow(china_prev),] %*% (W_late + 1e-9)
    assertthat::assert_that(all(dim(z[t_agg:nrow(z),])==dim(flight_prev_late)))
    llik[4] = sum(dpois(z[t_agg:nrow(z),], flight_prev_late, log=TRUE))

    if (isTRUE(visualise)) {
      par(mfrow=c(1,2))
      wuhan = cumsum(c(y_early[, 151], y_late[, 13]))
      beijing = cumsum(c(y_early[, 13], y_late[, 2]))
      plot(wuhan, main='Wuhan')
      lines(cumsum(exp_incr[, 151]), col=2)
      plot(beijing, main='Beijing')
      lines(cumsum(exp_incr[, 13]), col=2)
    }
    print(llik)
    sum(llik)
  }
}

