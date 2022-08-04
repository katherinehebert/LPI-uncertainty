# Script to make a little example plot(s) of accuracy and precision

library(ggplot2)
library(ggpubr)

set.seed(42)

theme_set(ggpubr::theme_transparent())

x = 1:10

df = data.frame(
  "x" = x,
  "y" = x*2
)
df$y_est = sapply(df$y, function(x) runif(1, x-5, x+5))

ggplot(df) +
  geom_line(aes(x = x, y = y), lty = 2) +
  coord_cartesian(ylim = c(0, 25))
ggsave("figures/presentation_demoAccuracy_truetrend.png", width = 6.41, height = 4.06)


ggplot(df) +
  geom_line(aes(x = x, y = y), lty = 2) +
  geom_smooth(aes(x = x, y = y_est), method = "lm", alpha = 0, col = "#303926", fill = "#cacfbc")+
  coord_cartesian(ylim = c(0, 25))
ggsave("figures/presentation_demoAccuracy_noci.png", width = 6.41, height = 4.06)


ggplot(df) +
  geom_line(aes(x = x, y = y), lty = 2) +
  geom_smooth(aes(x = x, y = y_est), method = "lm", col = "#303926", fill = "#cacfbc")+
  coord_cartesian(ylim = c(0, 25))
ggsave("figures/presentation_demoAccuracy.png", width = 6.41, height = 4.06)

png("figures/presentation_demoAccuracy_lambda.png")
x <- seq(-4, 4, length=100)
y <- dnorm(x)
plot(x,y, type = "l", lwd = 2, axes = FALSE, xlab = "", ylab = "")
abline(v = 0)
axis(1, at = 0, labels = c("λ"))
dev.off()

