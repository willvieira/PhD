#Plot "climate change" tendence from google scholar

#source
source("/Users/wvieira/Documents/GitHub/NativeFunctions/MyScatterPlot.R")

#data
sch <- read.table("/Users/wvieira/Documents/GitHub/PhD/litRev/data/scholar_dat.txt", head = T)
isi <- read.table("/Users/wvieira/Documents/GitHub/PhD/litRev/data/ISI_dat.txt", head = T)

#Plot
myplot(x = sch$year,
       y = log(sch$results),
       xlim =c(1940, 2015),
       ylim = c(1, 13.5),
       xlab = "time (year)",
       ylab = "frequency (log)",
       type = "l",
       newWindow = TRUE,
       width = 3.8,
       height = 3.3
)
lines(isi$year, log(isi$results), lty = 2)
legend("topleft", legend=c("Scholar", "ISI"), lty = 1: 2, cex = 0.95, bty = "n")

#export
library(extrafont)
loadfonts()
pdf("/Users/wvieira/Documents/GitHub/PhD/litRev/img/fig1.pdf", width = 3.8, height = 3.3, family="CM Roman")
myplot(x = sch$year,
       y = log(sch$results),
       xlim =c(1940, 2015),
       ylim = c(1, 13.5),
       xlab = "time (year)",
       ylab = "frequency (log)",
       type = "l",
       lwd = 2)
lines(isi$year, log(isi$results), lty = 2, lwd = 2)
legend("topleft", legend=c("Scholar", "ISI"), lty = 1: 2, cex = 0.95, bty = "n")
dev.off()
embed_fonts("/Users/wvieira/Documents/GitHub/PhD/litRev/img/fig1.pdf", outfile="/Users/wvieira/Documents/GitHub/PhD/litRev/img/fig1_em.pdf")
