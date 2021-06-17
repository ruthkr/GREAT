library(dplyr)
library(data.table)
library(ggplot2)

source("~/Downloads/Telegram Desktop/compare_registered_to_unregistered_model.R")
source("~/Downloads/Telegram Desktop/utils.R")
source("~/Downloads/Telegram Desktop/utils-calcs.R")
shifted.all.data.df <- readRDS("~/Downloads/Telegram Desktop/shifted.all.data.df.RDS")

plot_fit <- function(results_list, facets = c("arabidopsis", "brassica"), title = NULL) {
  data <- rbind(
    data.frame(
      shifted_time = results_list$ara.spline.data$shifted_time,
      y_truth = broom::augment(results_list$ara.fit)$mean.cpm,
      y_pred = broom::augment(results_list$ara.fit)$.fitted,
      type = facets[[1]]
    ),
    data.frame(
      shifted_time = results_list$bra.spline.data$shifted_time,
      y_truth = broom::augment(results_list$bra.fit)$mean.cpm,
      y_pred = broom::augment(results_list$bra.fit)$.fitted,
      type = facets[[2]]
    ),
    data.frame(
      shifted_time = results_list$combined.spline.data$shifted_time,
      y_truth = broom::augment(results_list$combined.fit)$mean.cpm,
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
      subtitle = paste0(
        "KS: ",
        "D=", round(results_list$ks$statistic, 2),
        ", p-value=", signif(results_list$ks$p.value, 3)
      ),
      x = "Shifted time",
      y = "mean.cpm"
    )
}

L_ara_ara <- compare_registered_to_unregistered_model(
  curr.sym = "BO4G120010",
  shifted.all.data.df,
  accessions = rep("Col0", 2)
)

L_bra_bra <- compare_registered_to_unregistered_model(
  curr.sym = "BO4G120010",
  shifted.all.data.df,
  accessions = rep("Ro18", 2)
)

L_ara_bra <- compare_registered_to_unregistered_model(
  curr.sym = "BO4G120010",
  shifted.all.data.df
)

list(
  plot_fit(L_ara_ara, title = "Ara/Ara - ", facets = c("arabidopsis1", "arabidopsis2")),
  plot_fit(L_bra_bra, title = "Bra/Bra - ", facets = c("brassica1", "brassica2")),
  plot_fit(L_ara_bra, title = "Ara/Bra - ")
) %>%
  patchwork::wrap_plots(ncol = 1, guides = "collect")
