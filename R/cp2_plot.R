# Loess smoothing model and contour plot code

library("dplyr")
library("ggplot2")
library("zoo")
library("metR")
library("viridis")


## Contour plot

list.files("data/output")

cp <- readRDS("data/output/cp_df.rds")

# fully processed

cp_gc <- filter(cp, sti == "gc")
cp_ct <- filter(cp, sti == "ct")

cp_gc <- filter(cp_gc, screening != Inf)
cp_ct <- filter(cp_ct, screening != Inf)

loess_gc <- loess(pia ~ or * screening, data = cp_gc, span = 0.25)
fit_gc <- expand.grid(list(or = seq(min(cp_gc$or), max(cp_gc$or), length.out = 100),
                           screening = seq(min(cp_gc$screening), max(cp_gc$screening), length.out = 100)))
fit_gc$pia <- as.numeric(predict(loess_gc, newdata = fit_gc))

loess_ct <- loess(pia ~ or * screening, data = cp_ct, span = 0.25)
fit_ct <- expand.grid(list(or = seq(min(cp_ct$or), max(cp_ct$or), length.out = 100),
                           screening = seq(min(cp_ct$screening), max(cp_ct$screening), length.out = 100)))
fit_ct$pia <- as.numeric(predict(loess_ct, newdata = fit_ct))

fit_gc$sti <- "NG"
fit_ct$sti <- "CT"

fit <- rbind(fit_gc, fit_ct)
fit <- rename(fit, PIA = pia)

ggplot(fit, aes(or, screening)) +
  geom_raster(aes(fill = PIA), interpolate = TRUE) +
  geom_contour(aes(z = PIA), col = "white", alpha = 0.5, lwd = 0.5) +
  theme_minimal() +
  theme(panel.spacing = unit(1.5, "lines")) +
  facet_grid(cols = vars(sti)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = "STI Screening Interval (Months)", x = "Per-Contact Medication Odds Ratio") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "B", direction = 1,
                     labels = scales::label_percent(1))
