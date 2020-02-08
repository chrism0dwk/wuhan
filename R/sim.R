#' Network SEIR metapopulation model
#'
#' @param N the population sizes
#' @param K_early China early between-cities connectivity
#' @param K_late China late between-cities connectivity
#' @param init_loc the name of the initial location
#' @param alpha the latent period
#' @param max_t the maximum time up to which to simulate
#' @return function which runs a simulation, with signature function(param).  Param is a vector of \code{c(beta, delta, gamma, I0W)}
#' @import deSolve
#' @export
NetworkODEModel = function(N, K_early, K_late, init_loc, alpha, max_t, t_control=23) {

  n = length(N) # number of cities

  ode_fn = function(t, y, parms) {

    ymat = matrix(y, nrow=n, byrow=FALSE)
    S = ymat[,1]
    E = ymat[,2]
    I = ymat[,3]
    R = ymat[,4]

    init_idx = which(rownames(K_early)==init_loc)
    K = K_early
    beta_wuhan = parms[1] # Early

    if (t >= t_control) {
      K = K_late
      beta_wuhan = parms[2]
    }

    # Within-city transmission vector
    beta = rep(parms[3], n)
    beta[init_idx] = beta_wuhan

    # Removal rate
    gamma = parms[4]

    lambda = beta * (I + (t(K/N)%*%I))
    dS = -lambda * S/N
    dE = lambda * S/N - alpha * E
    dI = alpha * E - gamma * I
    dR = gamma * I

    list(c(dS, dE, dI, dR))
  }

  t = 0:max_t

  simulate = function(param) {
    # Params: c(beta_w_early, beta_w_late, beta_w_china, gamma, I0W)
    I0W = param[5] # Initial infectives in Wuhan

    I0 = rep(0, length(N))
    I0[rownames(K_early)==init_loc] = I0W
    y = c(N-I0, rep(0, length(N)), I0, rep(0, length(N)))
    sim = deSolve::ode(y=y, times=t, func=ode_fn, parms=param)
    t=sim[,1]
    S=sim[,seq(2,len=n)]
    E=sim[,seq(2+n, len=n)]
    I=sim[,seq(2+2*n, len=n)]
    R=sim[,seq(2+3*n, len=n)]
    list(t=t, S=S, E=E, I=I, R=R)
  }
  simulate
}


NoiseGeneratingFunction = function(param, N, K, W, phi_mask=(rownames(K)=='Wuhan'),
                                   agg_up_to = 11, max_t=22, simulator=NA) {

  expected = simulator(param[1:3])
  p_detect = rep(1, length(N))
  p_detect[phi_mask] = param[4]

  # Increments on R
  exp_incr = t(t(diff(expected$R)) * p_detect)

  # Noise for China
  y_china = matrix(rpois(length(exp_incr), unlist(exp_incr)), ncol=ncol(exp_incr))

  # Noise for RoW
  china_prev = t(t(expected$I / N) * p_detect)
  flight_prev = china_prev %*% W
  y_row = matrix(rpois(length(flight_prev), unlist(flight_prev)), ncol=ncol(flight_prev))

  list(y=t(y_china), z=t(y_row))
}
