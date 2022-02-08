
## 04 Sim Min Exploratory Analysis

library("EpiModelHIV")
library("EpiModelHPC")

d <- merge_simfiles(simno = "3100", indir = "data/output/")
d

par(mar = c(3,3,1,1), mgp = c(2,1,0))
plot(d, y = "num", ylim = c(50000, 150000))
plot(d, y = "incid", ylim = c(0, 50))
plot(d, y = "ir100", ylim = c(0, 3))
plot(d, y = "i.num")
plot(d, y = "i.prev", ylim = c(0, 1))
plot(d, y = "ir100.sti")
plot(d, y = "new.aids.tot")

df <- as.data.frame(d)
head(df)
tail(df)

library(ggplot2)
ggplot() +
  geom_line(data = df, mapping = aes(time, i.prev, group = sim),
            alpha = 0.01, lwd = 0.25, color = "firebrick") +
  geom_bands(data = df, mapping = aes(time, i.num),
             lower = 0.1, upper = 0.9, fill = "firebrick") +
  ylim(c(0, 0.5)) +
  theme_minimal()
