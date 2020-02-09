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
                       W_early=W_late/31,
                       simulator=simulator)
  fit = optim(c(log(c(
    1.6, 1.6, 0.5, 0.5, 35
  )), logit(0.05)),
  llik,
  control = list(fnscale = -1),
  visualise = T)

  # Predictive simulation
  simulator_pred = NetworkODEModel(
    china_population$Population,
    K_early,
    K_late,
    'Wuhan',
    alpha = 1 / 4,
    max_t = 365 * 2,
    t_control = 23
  )
  sim = sim = simulator_pred(c(
    1.45520920, 1.65177535, 1.18430698, 0.57174427, 510.71167077, 0.01638659
  ))

  # Infection curves
  par(mfrow=c(1,2))
  china_agg = aggregate_by_province(sim$I, china_population$Province)
  plot(as.Date('2020-01-01')+sim$t[1:180], sim$I[1:180,151], type='l', main='Wuhan',
       xlab='Date', ylab=expression(I(t)))
  plot(as.Date('2020-01-01')+sim$t[1:180], rowSums(sim$I[1:180,]), main='China',
       xlab='Date', ylab=expression(I(t)), type='l')

  # Predictive plots
  plot_pred_detections = function(sim, horizon) {
    par(mfrow = c(2, 2))
    date_origin = as.Date('2020-01-01')
    china_prov = aggregate_by_province(t(china_cases), china_population$Province)
    pred_prov = aggregate_by_province(sim$R, china_population$Province)

    prov_names = c('Hubei', 'Beijing', 'Guangdong', 'Shanghai')
    for (prov in prov_names) {
      phi = ifelse(prov == 'Hubei', invlogit(fit$par[6]), 1)
      plot(
        date_origin + 1:35,
        cumsum(china_prov[, prov]),
        main = prov,
        pch = 20,
        xlab = 'Date',
        ylab = 'Cumulative detected',
        xlim = c(date_origin, date_origin + horizon),
        ylim = c(0, max(pred_prov[1:horizon, prov] * phi))
      )
      lines(date_origin + sim$t[1:horizon], pred_prov[1:horizon, prov] * phi, col =
              2)
    }
  }

  plot_pred_detections(sim, 40)
  plot_pred_detections(sim, 365)
}
