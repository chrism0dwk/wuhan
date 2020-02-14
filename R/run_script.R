#' Code run for analysis

run_script = function()
{
  # Fit model
  simulator = NetworkODEModel(china_population$Population,
                              K_early/31,
                              K_late/31,
                              'Wuhan',
                              alpha = 1 / 4,
                              35,
                              23)
  llik = LogLikelihood(y=china_cases,
                       z = world_cases,
                       pop_info=china_population,
                       K_early=K_early/31,
                       K_late=K_late/31,
                       W_early=W_early/31,
                       W_late=W_late/31,
                       sim_fun=simulator)
  fit = optim(c(log(c(
    1.6, 1.6, 0.5, 0.5, 35
  )), logit(0.05)),
  llik,
  control = list(fnscale = -1),
  visualise = T)

  params = c(exp(fit$par[1:5]), invlogit(fit$par[6]))


  # Predictive simulation
  simulator_pred = NetworkODEModel(
    china_population$Population,
    K_early/31,
    K_late/31,
    'Wuhan',
    alpha = 1 / 4,
    max_t = 365,
    t_control = 23
  )
  sim = simulator_pred(params)

  # Infection curves
  par(mfrow=c(1,2))
  china_agg = aggregate_by_province(sim$I, china_population$Province)
  plot(as.Date('2020-01-01')+sim$t[1:220], sim$I[1:220,151], type='l', main='Wuhan',
       xlab='Date', ylab=expression(I(t)))
  plot(as.Date('2020-01-01')+sim$t[1:220], rowSums(sim$I[1:220,]), main='China',
       xlab='Date', ylab=expression(I(t)), type='l')

  # Predictive plots
  plot_pred_detections = function(sim, horizon, log=FALSE) {
    par(mfrow = c(2, 2))
    date_origin = as.Date('2020-01-01')

    p_detect = rep(1, nrow(china_cases))
    p_detect[151] = invlogit(fit$par[6])
    expectedR = t(t(sim$R) * p_detect)

    china_prov = aggregate_by_province(t(china_cases), china_population$Province)
    pred_prov = aggregate_by_province(expectedR, china_population$Province)

    transform = identity
    ylab = 'Cumulative detected'
    if(isTRUE(log)) {
      transform = log10
      ylab = 'Log10 cumulative detected'
    }

    prov_names = c('Hubei', 'Beijing', 'Guangdong', 'Shanghai')
    for (prov in prov_names) {
      plot(
        date_origin + 1:35,
        transform(cumsum(china_prov[, prov])),
        main = prov,
        pch = 20,
        xlab = 'Date',
        ylab = ylab,
        xlim = c(date_origin, date_origin + horizon),
        ylim = c(0, max(transform(pred_prov[1:horizon, prov])))
      )
      lines(date_origin + sim$t[1:horizon], transform(pred_prov[1:horizon, prov]), col=2)
    }
  }

  plot_pred_detections(sim, 40)
  plot_pred_detections(sim, 365, log=T)


  # Noise function for resampling
  noise = NoiseGeneratingFunction(params, china_population$Population,
                                  W_early, W_late)
  llik1 = LogLikelihood(
    y = noise$y,
    z = noise$z,
    pop_info = china_population,
    K_early = K_early / 31,
    K_late = K_late / 31,
    W_early = W_early / 31,
    W_late = W_late / 31,
    sim_fun = simulator
  )
  fit_check = optim(fit$par,
                    llik1,
                    control = list(fnscale = -1),
                    visualise = T)
}


extract_peak_time = function(sim, city) {
  sim$t[which.max(sim$I[,city])]
}

plot_I_curves = function(sims) {

  plot(as.Date('2020-01-01'),0, xlim=as.Date('2020-01-01')+c(0, 160), ylim=c(0, 8.5e5), type='n', ylab=expression(I(t)), xlab='Date', main='Wuhan')
  sapply(sims, function(sim) lines(as.Date('2020-01-01')+sim$t, sim$I[,151], col=rgb(1, 0, 0, 0.1)))

  plot(as.Date('2020-01-01'),0, xlim=as.Date('2020-01-01')+c(0, 160), ylim=c(0, 8.5e5), type='n', ylab=expression(I(t)), xlab='Date', main='Wuhan')
  sapply(sims, function(sim) lines(as.Date('2020-01-01')+sim$t, sim$I[,151], col=rgb(1, 0, 0, 0.1)))
}

province_infec_curve = function(sims, prov='Hubei') {
  w = as.data.frame(lapply(sims, function(sim) {
    aggE = aggregate_by_province(sim$E, china_population$Province)
    aggI = aggregate_by_province(sim$I, china_population$Province)
    aggE[, prov] + aggI[, prov]
  }))
  t(apply(w, 1, function(x) c(mean(x), quantile(x, probs=c(0.025, 0.975)))))
}

china_infec_curve = function(sims) {
  china = as.data.frame(lapply(sims, function(sim) rowSums(sim$I)+rowSums(sim$E)))
  t(apply(china, 1, function(x) c(mean(x), quantile(x, probs=c(0.025, 0.975)))))
}


uk_predicted_imports = function(sim) {
  china_prev = t(t(sim$I)/china_population$Population)

  early = china_prev[1:22,] %*% (W_early[,'United.Kingdom']/31)
  late = china_prev[23:nrow(china_prev),] %*% (W_late[,'United.Kingdom']/31)
  uk_imports = rbind(early, late)
}

uk_imports_l = function(sims) {
  sapply(sims, uk_predicted_imports)
}
