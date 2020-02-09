#' Code run for analysis

run_script = function()
{

  # Fit model
  simulator = NetworkODEModel(china_population$Population,
                              K_early, K_late, 'Wuhan', alpha=1/4, 35, 23)
  llik = LogLikelihood(china_cases, z=world_cases, china_population,
                       K_early, K_late, W_early, W_late, simulator)
  fit = optim(c(log(c(1.6, 1.6, 0.5, 0.5, 35)), logit(0.05)),
              llik, control=list(fnscale=-1), visualise=T)

  # Predictive simulation
  simulator_pred = NetworkODEModel(china_population$Population, K_early, K_late,
                                   'Wuhan', alpha=1/4, max_t=365*2, t_control=23)
  sim = sim = simulator_pred(c(4.87856467, 2.05144661, 0.93721705, 0.91207631, 0.06335245, 0.02623697))


  # Predictive plots
  plot_pred_detections = function(sim, horizon) {
    par(mfrow=c(2,2))
    date_origin = as.Date('2020-01-01')
    china_prov = aggregate_by_province(t(china_cases), china_population$Province)
    pred_prov = aggregate_by_province(sim$R, china_population$Province)

    prov_names = c('Hubei', 'Beijing', 'Guangdong', 'Shanghai')
    for (prov in prov_names) {
      phi = ifelse(prov == 'Hubei', 0.026, 1)
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
