plot_fit <- function(results_list, facets = c("arabidopsis", "brassica"), title = NULL) {
  data <- rbind(
    data.frame(
      shifted_time = results_list$ara.spline.data$shifted_time,
      y_truth = broom::augment(results_list$ara.fit)$mean_cpm,
      y_pred = broom::augment(results_list$ara.fit)$.fitted,
      type = facets[[1]]
    ),
    data.frame(
      shifted_time = results_list$bra.spline.data$shifted_time,
      y_truth = broom::augment(results_list$bra.fit)$mean_cpm,
      y_pred = broom::augment(results_list$bra.fit)$.fitted,
      type = facets[[2]]
    ),
    data.frame(
      shifted_time = results_list$combined.spline.data$shifted_time,
      y_truth = broom::augment(results_list$combined.fit)$mean_cpm,
      y_pred = broom::augment(results_list$combined.fit)$.fitted,
      type = "combined"
    )
  )

  data %>%
    tidyr::pivot_longer(
      cols = -c(shifted_time, type),
      names_to = "y"
    ) %>%
    ggplot() +
    aes(x = shifted_time, y = value, color = y) +
    geom_point(alpha = 0.5, shape = 1) +
    geom_line() +
    facet_wrap(~type, nrow = 1) +
    labs(
      title = paste0(
        title,
        "AIC=",
        round(results_list$seperate.AIC), ":", round(results_list$combined.AIC),
        ", BIC=",
        round(results_list$seperate.BIC), ":", round(results_list$combined.BIC),
        " (sep:reg)"
      ),
      # subtitle = paste0(
      #   "KS: ",
      #   "D=", round(results_list$ks$statistic, 2),
      #   ", p-value=", signif(results_list$ks$p.value, 3)
      # ),
      x = "Shifted time",
      y = "mean_cpm"
    ) +
    scale_x_continuous(limits = c(10, 60))
}
