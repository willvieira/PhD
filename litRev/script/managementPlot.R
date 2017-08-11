#Hypothesized interaction between forest management and ecological processes

#packages
library(pBrackets)
library(extrafont)
loadfonts()

#generate management data
fm <- seq(0, 10, 0.05)

#models
linFunc <- function(x, b) x*b + 2

#plot
pdf("/Users/wvieira/GitHub/PhD/litRev/img/fig2.pdf", width = 8.5, height = 2.2, family = "serif")
par(mfrow = c(1, 4), mar = c(2.25, 1.7, .5, 1.5))

#Plantation
b = 2.7
plot(fm, linFunc(fm, b), ylim = c(-4, 32), type = "l", lwd = 3, xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty='l')
mtext(side = 3, expression(alpha), line = -2, at = 10.7)
abline(1.8,0, lty = 3)
brackets(10.2, -5.5, 10.2, 1.8, ticks = 0.5, curvature = 0.5, type = 1, lwd = 1, lty = 1)
mtext(side = 1, "Natural base process", line = - 1.9, at = 4.8, cex = 0.9)
mtext(side = 1, "Plantation", line = 0.8, cex = 1)
#Thinning
b = 2
plot(fm, linFunc(fm, b), ylim = c(-4, 32), type = "l", lwd = 3, xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty='l')
mtext(side = 3, expression(theta), line = -4.7, at = 10.7)
abline(1.8,0, lty = 3)
mtext(side = 1, "Thinning", line = 0.8, cex = 1)
#Enrichment
b = 1.25
plot(fm, linFunc(fm, b), ylim = c(-4, 32), type = "l", lwd = 3, xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty='l')
mtext(side = 3, expression(beta), line = -7.4, at = 10.7)
abline(1.8,0, lty = 3)
mtext(side = 1, "Enrichment planting", line = 0.8, cex = 1)
#Disturbance
b = 3.1
plot(fm, linFunc(fm, b), ylim = c(-4, 32), type = "l", lwd = 3, xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty='l')
mtext(side = 3, expression(epsilon), line = - 0.5, at = 10.7)
abline(1.8,0, lty = 3)
mtext(side = 1, "Harvest", line = 0.8, cex = 1)

#Close window
dev.off()

#Embed fonts into the file
embed_fonts("/Users/wvieira/GitHub/PhD/litRev/img/fig2.pdf", outfile="/Users/wvieira/GitHub/PhD/litRev/img/fig2_em.pdf")
