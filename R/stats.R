#' Bootstrap for confidence intervals
#'
#' @param param_hat parameters c(beta, gamma, I0W, phi) on the log(it) scale
#' @param K the flight connectivity matrix within China
#' @param W the flight connectivity matrix outside China
#' @param alpha the latent period
#' @param max_t time to simulate up to
#' @return a data.frame of bootstrapped parameter samples
#' @export
bootstrap = function(param_hat, K, W, alpha, max_t, n_samples=10) {

  simulator = NetworkODEModel(china_population, K, 'Wuhan', alpha, max_t)

  bs = function() {
    sim = NoiseGeneratingFunction(param=param_hat, N=china_population, K=K, W=W, simulator=simulator)

    llik = LogLikelihood(y=sim$y, z=sim$z, N=china_population, K=K, W=W, sim_fun=simulator)
    p_hat = optim(log(param_hat), llik, control=list(fnscale=-1, maxit=2000), visualise=TRUE)
    p_hat
  }

  n_sim = 0
  samples = NULL
  while(n_sim < n_samples) {
    fit = bs()
    if(fit$convergence == 0) {
      samples = rbind(samples, fit$par)
      n_sim = n_sim + 1
    }
  }
  samples = as.data.frame(samples)
  names(samples) = c('beta','gamma','I0W','phi')
  samples
}

#' Calculates 95% CI around median of an empirical distn.
#' @export
quantile_ci = function(x) quantile(x, probs=c(0.025, 0.5, 0.975))


#' Summarises samples from bootstrapping
#'
#' @param samples an \eqn{n \times 4} \code{data.frame} of \code{beta}, \code{gamma}, \code{I0W}, \code{phi}
#' @return table of summary data
#' @export
summary_samples = function(samples) {
  beta = exp(samples$beta)
  gamma = exp(samples$gamma)
  I0W = exp(samples$I0W)
  phi = invlogit(samples$phi)

  res = rbind(quantile_ci(beta),
              quantile_ci(gamma),
              quantile_ci(1/gamma),
              quantile_ci(I0W),
              quantile_ci(phi),
              quantile_ci(beta/gamma))
  rownames(res) = c('beta','gamma','lat_period','I0W','phi','R0')
  colnames(res) = c('2.5%','50%','97.5%')
  as.data.frame(res)
}
