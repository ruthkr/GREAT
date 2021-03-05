# Functions ----

add_n_points <- function(vect, n) {
  a <- vect[[1]]
  b <- vect[[2]]

  divisions <- seq(0, 1, 1 / (n + 1))
  distance <- b - a
  new_vect <- a + divisions * distance

  return(c(a, new_vect[2:(n + 1)]))
}

expand_x_space <- function(data, x_var, n) {
  x_vect <- data %>%
    dplyr::pull({{ x_var }}) %>%
    sort()

  # Matrix of consecutive pairs
  x_pairs <- cbind(x_vect[-length(x_vect)], x_vect[-1])

  # Expand space
  new_x <- x_pairs %>%
    split(., row(.)) %>%
    purrr::map(
      function(x) {
        return(add_n_points(x, n))
      }
    ) %>%
    purrr::reduce(c)

  return(c(new_x, x_vect[[length(x_vect)]]))
}

expand_y_space <- function(data, x_var, y_var, n) {

  new_data <- data.frame(
  x = expand_x_space(data, {{ x_var }}, n)
  ) %>%
  dplyr::rename({{x_var}} := x)

  new_data <- new_data %>%
    dplyr::left_join(
      data,
      by = quo_name(enquo(x_var))
    ) %>%
    imputeTS::na_interpolation(option = "linear")

  return(new_data)
}

# Test ----

# set.seed(12)
#
# data <- data.frame(
#   x = c(1, 2, 3, 4, 6, 7),
#   y = runif(6)
# )
#
# data %>% plot()
# expand_y_space(data, x, y, 10) %>% plot()
#
#
# # Real data ----
#
# shifted.all.data.df <- readRDS("~/Downloads/Telegram Desktop/shifted.all.data.df.RDS")
#
# new_data <- shifted.all.data.df %>%
#   dplyr::filter(accession == "Col0") %>%
#   dplyr::select(shifted.time, mean.cpm) %>%
#   expand_y_space(shifted.time, mean.cpm, 0) %>%
#   distinct()
#
#
# num.spline.params <- 6 # number of parameters for each spline fitting (degree and this used to calculate num knots).
# # num.registration.params <- 2 # stretch, shift
# # num.obs <- nrow(combined.spline.data)
# stats::lm(
#   mean.cpm ~ splines::bs(shifted.time, df = num.spline.params, degree = 3),
#   # data = ara.spline.data
#   data = new_data
#   # data = shifted.all.data.df %>% dplyr::filter(accession == "Col0")
# ) %>%
#   BIC()
