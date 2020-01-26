#' Network SEIR metapopulation model
#'
#' @param N the population sizes
#' @param K a contact matrix
#' @param init_loc the name of the initial location
#' @param alpha the latent period
#' @return list of matrices containing simulated results
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
    sim = ode(y=y, times=t, func=func, parms=param)
    t=sim[,1]
    S=sim[,seq(2,len=n)]
    E=sim[,seq(2+n, len=n)]
    I=sim[,seq(2+2*n, len=n)]
    R=sim[,seq(2+3*n, len=n)]
    list(t=t, S=S, E=E, I=I, R=R)
  }

}
