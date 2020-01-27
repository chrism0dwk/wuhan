#' Network SEIR metapopulation model
#'
#' @param N the population sizes
#' @param K a contact matrix
#' @param init_loc the name of the initial location
#' @param alpha the latent period
#' @return list of matrices containing simulated results
#' @import deSolve
#' @export
NetworkODEModel = function(N, K, init_loc, alpha, max_t) {

  n = length(N) # number of cities

  func = function(t, y, parms) {
    ymat = matrix(y, nrow=n, byrow=FALSE)
    S = ymat[,1]
    E = ymat[,2]
    I = ymat[,3]
    R = ymat[,4]
    beta = parms[1]
    gamma = parms[2]

    lambda = beta * (I/N + (t(K/N)%*%I/N))
    dS = -lambda * S
    dE = lambda * S - alpha * E
    dI = alpha * E - gamma * I
    dR = gamma * I

    list(c(dS, dE, dI, dR))
  }

  t = 0:max_t


  simulate = function(param) {
    I0W = param[3] # Initial infectives in Wuhan
    I0 = rep(0, length(N))
    I0[rownames(K)==init_loc] = I0W
    y = c(N-I0, rep(0, length(N)), I0, rep(0, length(N)))
    sim = deSolve::ode(y=y, times=t, func=func, parms=param)
    t=sim[,1]
    S=sim[,seq(2,len=n)]
    E=sim[,seq(2+n, len=n)]
    I=sim[,seq(2+2*n, len=n)]
    R=sim[,seq(2+3*n, len=n)]
    list(t=t, S=S, E=E, I=I, R=R)
  }

}


<<<<<<< HEAD
NoiseGeneratingFunction = function(param, N, K, W, phi_mask=rownames(K)=='Wuhan',
                                   agg_up_to = 11, max_t=22, r=2, simulator=None) {
=======
NoiseGeneratingFunction = function(param, N, K, W, phi_mask=(rownames(K)=='Wuhan'),
                                   agg_up_to = 11, max_t=22, simulator=None) {
>>>>>>> master

  expected = simulator(param[1:3])
  p_detect = rep(1, length(N))
  p_detect[phi_mask] = param[4]

  # Increments on R
  exp_incr = t(t(diff(expected$R)) * p_detect)

  # Noise for China
<<<<<<< HEAD
  y_china = matrix(rnbinom(length(exp_incr), size=r, mu=unlist(exp_incr)), ncol=ncol(exp_incr))
=======
  y_china = matrix(rpois(length(exp_incr), unlist(exp_incr)), ncol=ncol(exp_incr))
>>>>>>> master

  # Noise for RoW
  china_prev = t(t(expected$I / N) * p_detect)
  flight_prev = china_prev %*% W
<<<<<<< HEAD
  y_row = matrix(rnbinom(length(flight_prev), size=r, mu=unlist(flight_prev)), ncol=ncol(flight_prev))
=======
  y_row = matrix(rpois(length(flight_prev), unlist(flight_prev)), ncol=ncol(flight_prev))
>>>>>>> master

  list(y=t(y_china), z=t(y_row))
}
