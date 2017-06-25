#Hypothesized interaction between forest management and ecological processes

#packages
library(pBrackets)
library(extrafont)
loadfonts()

#generate management data
fm <- seq(0, 10, 0.05)

#models
alpha <- function(x) x^(2/6) + 2
beta  <- function(x) x*3 + 2
theta <- function(x) x*1.5 + 2
eps   <- function(x) x^(2/6) + 2

#plot
pdf("/Users/wvieira/GitHub/PhD/litRev/img/fig2.pdf", width = 7.97, height = 2.6, family = "serif")
par(mfrow = c(1, 3), mar = c(2.25, 1.7, .5, 1.5))

#Plantation
plot(fm, alpha(fm), ylim = c(min(alpha(fm)) - 0.4, max(alpha(fm))), type = "l", lwd = 3, xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty='l')
mtext(side = 3, expression(alpha), line = -1, at = 10.7)
abline(1.99,0, lty = 3)
brackets(10.2, 1.51, 10.2, 1.98, ticks = 0.5, curvature = 0.5, type = 1, lwd = 1, lty = 1)
mtext(side = 1, "Natural base process", line = - 2.05, at = 6.3, cex = 1)
mtext(side = 1, "Plantation", line = 0.8, cex = 1)
#Thinning
plot(fm, beta(fm), ylim = c(min(beta(fm)) - 6, max(beta(fm))), type = "l", lwd = 3, xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty='l')
points(fm, theta(fm), ylim = c(min(theta(fm)) - 6, max(theta(fm))), type = "l", lwd = 3, xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty='l')
mtext(side = 3, expression(theta), line = -1, at = 10.7)
mtext(side = 3, expression(beta), line = -8, at = 10.7)
abline(1.8,0, lty = 3)
mtext(side = 1, "Thinning", line = 0.8, cex = 1)
#Disturbance
plot(fm, eps(fm), ylim = c(min(eps(fm)) - 0.4, max(eps(fm))), type = "l", lwd = 3, xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty='l')
mtext(side = 3, expression(epsilon), line = -1, at = 10.7)
abline(1.99,0, lty = 3)
mtext(side = 1, "Harvest", line = 0.8, cex = 1)

#Close window
dev.off()

  #Embed fonts into the file
embed_fonts("/Users/wvieira/GitHub/PhD/litRev/img/fig2.pdf", outfile="/Users/wvieira/GitHub/PhD/litRev/img/fig2_em.pdf")
