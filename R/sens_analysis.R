#' Sensitivity analysis for alpha

par_init_ = c(0.4, 0.1428571, 1.0, 0.5)

#' Fits the nCoV model
#'
#' @param alpha the assumed latent period
#' @param K the flight connectivity matrix within China
#' @param W the flight connectivity matrix internationally
#' @param init_loc the initial city within China
#' @param max_t the time up to which to make inference
#' @return output from optim()
#' @export
fit_model = function(alpha=1/4, K=NA, W=NA, init_loc='Wuhan', max_t=22) {

  simulator = NetworkODEModel(china_population, K, 'Wuhan', alpha, max_t)
  llik = LogLikelihood(china_cases[,1:max_t], world_cases[,1:max_t],china_population, K, W, simulator)
  p_hat = optim(log(par_init_), llik, control=list(fnscale=-1, maxit=2000), visualise=TRUE)
  p_hat
}


#' Sensitivity to alpha
#'
#' @param alpha a vector of alpha parameters to explore sensitivity to
#' @param K the flight connectivity matrix within China
#' @param W the flight connectivity matrix internationally
#' @return a list of optim outputs and alphas
#' @export
alpha_sensitivity = function(alpha=1/c(3.6, 4.0, 4.4, 5.0, 6.0), K, W) {
  res = list()
  for (i in 1:length(alpha)) {
    res[[i]] = fit_model(alpha=alpha[i], K=K, W=W)
  }
  res$alpha=alpha
  res
}
