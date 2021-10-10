## A noisy sine wave as query
library(dtw)
idx<-seq(0,6.28,len=100);
query<-sin(idx)+runif(100)/10;

## A cosine is for template; sin and cos are offset by 25 samples
# idx_2 <- idx/2
template<- cos(idx)

## Find the best match with the canonical recursion formula
library(dtw);
alignment<-dtw(query,template,keep=TRUE);

## Display the warping curve, i.e. the alignment curve
plot(alignment,type="threeway")

## Align and plot with the Rabiner-Juang type VI-c unsmoothed recursion
plot(
  dtw(query,template,keep=TRUE,
      step=rabinerJuangStepPattern(6,"c")),
  type="twoway",offset=-2);

# Plotting two-ways
dtwPlotTwoWay(alignment)

## See the recursion relation, as formula and diagram
rabinerJuangStepPattern(6,"c")
plot(rabinerJuangStepPattern(6,"c"))

library(dplyr)

data_test <- data.frame(
  idx = idx,
  query = query,
  template = template
) %>%
  tidyr::pivot_longer(cols = c("query", "template"),
                      names_to = "group")


ggplot2::ggplot(data = data_test) +
  ggplot2::aes(x = idx, y = value, colour = group) +
  ggplot2::geom_line() +
  ggplot2::geom_point()
