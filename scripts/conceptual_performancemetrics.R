# figures for conceptual diagram of measurements

library(dplyr)
library(ggplot2)

theme_set(ggpubr::theme_pubr())


## Accuracy:

x <- seq(from = 1, to = 1.5, length.out = 10)
v <- seq(from = 1, to = 1.1, length.out = 10)

df <- data.frame(
  "time" = 1:10,
  "LPI_final" = x,
  "LPI_final_true" = v
)
df$diff <- df$LPI_final - df$LPI_final_true

ggplot(df, aes(x = time)) +
  geom_line(aes(y = LPI_final), col = "black", lwd = 1) +
  geom_line(aes(y = LPI_final_true), lty = 2, col = "black", lwd =1) +
  ggforce::geom_link(data = df[-1,],
                     aes(x = time, y = LPI_final, 
                         xend = time, yend = LPI_final_true, 
                         colour = diff),
                     arrow = arrow(ends = "both", type = "closed"), size = 1) +
  labs(x = "", y = "Living Planet Index")  +
  scale_color_viridis_c(option = "plasma", end = .8)

## Precision: success rate

x <- seq(from = 1, to = 1.5, length.out = 10)
v <- c(seq(from = 1, to = 1.1, length.out = 3),
       seq(from = 1.15, to = 1.35, length.out = 7))

df <- data.frame(
  "time" = 1:10,
  "LPI_final" = x,
  "LPI_final_true" = v,
  "CI_low" = x - seq(0, 0.1, length.out = 10),
  "CI_high" = x + seq(0, 0.1, length.out = 10)
)
df$Success = c(rep("Success", 5), rep("Failure", 5))

ggplot(df, aes(x = time)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = .4, fill = "grey") +
  geom_line(aes(y = LPI_final), col = "black", lwd = 1) +
  geom_line(aes(y = LPI_final_true), lty = 2, col = "black", lwd = 1) +
  geom_point(aes(y = LPI_final_true, col = Success), size = 4) +
  # ggforce::geom_link(data = df[-1,],
  #                    aes(x = time, y = LPI_final, 
  #                        xend = time, yend = LPI_final_true, 
  #                        colour = diff),
  #                    arrow = arrow(ends = "both", type = "closed"), size = 1) +
  labs(x = "", y = "Living Planet Index")  +
  scale_color_manual(values = c("red", "#4bc96c"), labels = c("Failure", "Success")) +
  theme(legend.position = "right")

## Precision: CI diff

x <- seq(from = 1, to = 1.5, length.out = 10)
v <- seq(from = 1, to = 1.45, length.out = 10)
z <- seq(from = 1, to = 1.4, length.out = 10)

df <- data.frame(
  "time" = 1:10,
  "LPI_final" = x,
  "CI_low" = x - seq(0, 0.1, length.out = 10),
  "CI_high" = x + seq(0, 0.1, length.out = 10),
  "LPI_final_true" = v,
  "LPI_raw" = x,
  "CI_low_raw" = x - seq(0, 0.13, length.out = 10),
  "CI_high_raw" = x + seq(0, 0.13, length.out = 10)
)

ggplot(df, aes(x = time)) +
  geom_ribbon(aes(ymin = CI_low_raw, ymax = CI_high_raw), alpha = .1, fill = "red") +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 1, fill = "grey86") +
  geom_line(aes(y = LPI_final), col = "black", lwd = .5) +
  geom_ribbon(aes(ymin = CI_low_raw, ymax = CI_high_raw), alpha = 0, fill = "red", lwd = .5, lty = 2, col = "red") +
#  geom_line(aes(y = LPI_final_true), lty = 2, col = "black", lwd = .5) +
  ggforce::geom_link(aes(x = time, y = CI_low, xend = time, yend = CI_low_raw), col = "blue", lwd = 1) +
  ggforce::geom_link(aes(x = time, y = CI_high, xend = time, yend = CI_high_raw), col = "blue", lwd = 1) +
  labs(x = "", y = "Living Planet Index") 

x <- seq(from = 1, to = 1.5, length.out = 10)
v <- seq(from = 1, to = 1.45, length.out = 10)

df <- data.frame(
  x = rnorm(100, 1.3, 0.1),
  y = rnorm(100, 1.5, 0.1)
)

ggplot(df) +
  geom_density(aes(x = x), col = "grey70", fill = "grey", alpha =.1, lty = 2) +
  geom_density(aes(x = y), col = "grey70", fill = "grey", alpha =.1) +
  geom_vline(aes(xintercept = quantile(df$x, 0.025)), col = "blue", lty = 2) +
  geom_vline(aes(xintercept = quantile(df$x, 0.975)), col = "red", lty = 2) +
  geom_vline(aes(xintercept = quantile(df$y, 0.025)), col = "blue") +
  geom_vline(aes(xintercept = quantile(df$y, 0.975)), col = "red") +
  geom_segment(aes(x = quantile(df$x, 0.025), xend = quantile(df$y, 0.025),
                   y = 4, yend = 4), arrow = arrow(ends = "both"), col = "blue") +
  geom_segment(aes(x = quantile(df$x, 0.975), xend = quantile(df$y, 0.975),
                   y = 4, yend = 4), arrow = arrow(ends = "both"), col = "red")
  
