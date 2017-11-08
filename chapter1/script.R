#'---
#'title: "Linking State Transition Models and forest management"
#'author: Willian Vieira
#'date: "June, 28 2017 - last update: September, 17 2017"
#'output:
#'    html_document:
#'      toc: true
#'---
#'<!--Rscript -e "rmarkdown::render('script.R')" -->
#'
#'# Introduction
#' Here I am using a four state transition model (STM) to understand how forest management
#' can increase forest resilience.
#' The model has four process-parameters that define the transition rate between the
#' four states: Regeneration, Boreal, Temperate and Mixed (Figure 1).
#'
#' <center><img src="img/fig1.png" alt="Figure 1" style="width: 550px;"/></center>
#'
#' Each parameter represents a different process: $\alpha$ is the colonisation
#' rate (R -> B, T, M), $\beta$ the succession (B, T -> M), $\theta$ the exclusion
#' extintion (M -> B, T) and $\epsilon$ the pertubation (B, T, M -> R).
#'
#'# Example of running model
#'
# call the model
source("vissault_model_v3.R")
#'
#'#################################
#'## Running the model to equilibrium {#ov}
#'#################################

# parameters
params = read.table("pars.txt", row.names = 1)

# Wrapper to collect parameters for a given set of environmental conditions
pars = get_pars(ENV1 = 0, ENV2 = 0, params, int = 3)
pars
#'
# Get the transition matrix
transition_matrix = get_matrix(ENV1 = 0, ENV2 = 0, int = 3)
transition_matrix$MAT
#'
# Equilibrium of the model for the specific parameters (pars)
eq = get_eq(pars)
eq
#' `$eq` is the proportion of each state at equilibrium and `$ev` the eigenvalue.
#'
#'#################################
#'## See the model behavior graphically
#'#################################
#'
#' Starting with 25% of the proportion for each state, we can graphically see the
#' behavior of the model until reaching equilibrium.

# data frame
time <- seq(1, 50, by = 0.1)
dm2 <- data.frame(matrix(NA, nrow = length(time), ncol = 4))

# Initial condition
y = c(B = 0.25, T = 0.25, M = 0.25)
dm2[1, ] = c(y, 0.25)
#'
# Function to get behavior
behavior <- function(ENV1, ENV2) {

  # paramters depending on environment
  pars = get_pars(ENV1 = ENV1, ENV2 = ENV2, params, int = 3)

  # getting eq
  for(i in 2:length(time)) {
    eq = runsteady(y = y, func = model, parms = pars, times = c(0, i))[[1]]
    eq[4] <- 1 - eq[1] - eq[2] - eq[3]
    dm2[i, ] <- eq
  }

  # plot
  par(family = 'serif', cex = 0.8, mar = c(4, 4, 1, 1))
  plot(time, dm2[, 1],
       type = "l",
       lwd = 1.5,
       ylim = c(0, 1),
       ylab = "",
       xlab = "")
  lines(dm2[, 2], col = 2, lwd = 1.5)
  lines(dm2[, 3], col = 3, lwd = 1.5)
  lines(dm2[, 4], col = 4, lwd = 1.5)
}
#'
#+ fig.width = 9.5, fig.height = 8
Ev1 <- c(-1.75, -1.5, -1.25, -1, -0.75, -0.5, 0, 0.5, 1)

par(mfrow = c(3, 3))
for(i in 1: length(Ev1)) {
  behavior(ENV1 = Ev1[i], ENV2 = 0)
  mtext(paste("ENV1 = ", Ev1[i], sep = ""), 3, cex = 0.85)
  abline(h = seq(0, 1, 0.2), lty = 3, col = "gray")
  if(i == 4) mtext("State proportion", 2, line = 2, cex = 0.8)
  if(i == 8) mtext("Time (year*5)", 1, line = 2, cex = 0.8)
  if(i == 1) legend(35, 0.45, c("B", "T", "M", "R"), col = 1:4, lty = 1, bty = "n", cex = 0.8)
}
#'
#'##################################
#'# Model behavior depending on climate
#'#################################
#'
#' Now, let's take a look in the parameter variation depending on the climate envelop.
#'
# Varing temperature (ENV1) from -2 to 2 and then plot each parameter variation
tempVar <- seq(-2, 2, length.out = 100)
dat.Par <- data.frame(matrix(NA, ncol = length(pars) - 2, nrow = length(tempVar)))
names(dat.Par) <- names(pars)[1:7]

#loop
for(i in 1: length(tempVar)) {
  dat.Par[i, ] = get_pars(ENV1 = tempVar[i], ENV2 = 0, params, int = 3)[1:7]
}

# plot
#+ fig.width = 6, fig.height = 5.5
plot(tempVar, dat.Par[, 1], type = "l", lwd = 1.5,
     ylab = "parameter value", xlab = "Temperature variation", ylim = c(0, 1))
for(i in 2: 7) {
  points(tempVar, dat.Par[, i], type = "l", col = i, lwd = 1.5)
}
legend("topright", names(pars)[1:7], col = 1:7, lwd = 1.5, bty = "n")
abline(v = c(-1.35, -.65, 0.45), lty = 3, col = "gray")
mtext(c("Boreal", "Mixed", "Temperate"), 3, at = c(-1.35, -0.65, 0.45))
#'
#' With an idea of the parameter variation depending on the climate, we can see the model equilibrium
#' depending on the climate.
#'
# Varing temperature (ENV1) from -2 to 2 and then plot the eq proportion + eigenvalue
dat.eq <- data.frame(matrix(NA, ncol = 6, nrow = 100))
names(dat.eq) <- c("temp", "B", "T", "M", "R", "ev")
dat.eq$temp <- seq(-2, 1.2, length.out = 100)

for(i in 1: dim(dat.eq)[1]) {
  pars = get_pars(ENV1 = dat.eq$temp[i], ENV2 = 0, params, int = 3)
  # option to add management (change 0 for a value of management intensity)
  pars[c(2, 6, 7)] <- c((pars[2]*0) + pars[2], (pars[6]*0) + pars[6], (pars[7]*0) + pars[7])
  eq <- get_eq(pars)
  dat.eq[i, 2] <- eq$eq[1]
  dat.eq[i, 3] <- eq$eq[2]
  dat.eq[i, 4] <- eq$eq[3]
  dat.eq[i, 5] <- 1 - eq$eq[1] - eq$eq[2] - eq$eq[3]
  dat.eq[i, 6] <- eq$ev
}

# plot
#+ fig.width = 7, fig.height = 5.5
par(mar = c(4.5, 4.5, 2, 4.5))
plot(dat.eq$temp, dat.eq$B, type = "l", lwd = 1.5, ylim = c(0, 1), xlab = "Latitudinal gradient", ylab = "Occupancy")
for(i in 3: 5) points(dat.eq$temp, dat.eq[, i], type = "l", lwd = 1.5, col = i - 1)
par(new = T)
plot(dat.eq$temp, dat.eq$ev, type = "l", lty = 2, axes = F, xlab = NA, ylab = NA)
axis(side = 4)
mtext(side = 4, line = 3, 'Largest real part')
legend(-.8, 0, legend = names(dat.eq[-1]), lwd = 1.5, col = 1:4, lty = c(1, 1, 1, 1, 2), bty = "n")
abline(v = c(-1.54, -1.09, 0.24, 0.53), lty = 3, col = "gray")
mtext(c("Boreal", "Mixed", "Temperate"), 3, at = c(-1.35, -0.65, 0.45))

#'#################################
#'# Add hypothetical disturbance
#'#################################
#'
#' Setting up the data frame with two different environmental conditions is a way to
#' see the behavior of the model after a disturbance.
#' In this case, the model will start with 25% of the proportion for each state and after
#' reaching equilibrium, the environmnetal patters will change.
#'
# Function to produce the data frame with the equilibrium proportion of each state
# based on the environmental variation.
behavior <- function(envComb1, envComb2) {
	xx <- seq(1, 40, by = 1)
	dat <- data.frame(matrix(NA, nrow = length(xx), ncol = 6))
	names(dat) <- c("ENV1", "ENV2", "B", "T", "M", "R")
	dat[c(1:20), c(1,2)] <- envComb1 #environmental 1
	dat[c(21: dim(dat)[1]), c(1,2)] <- envComb2 #environmental 2

	# Initial condition
	dat[1, c(3:6)] = c(B = 0.25, T = 0.25, M = 0.25, R = 0.25)

	# Behavior
	for(i in 2:length(xx)) {
		pars = get_pars(ENV1 = dat[i, 1],ENV2 = dat[i, 2], params, int = 3)
		y = c(B = dat[(i - 1), 3], T = dat[(i - 1), 4], M = dat[(i - 1), 5])
		eq = runsteady(y = y, func = model, parms = pars, times = c(0, i))[[1]]
	  eq[4] <- 1 - eq[1] - eq[2] - eq[3]
	  dat[i, c(3:6)] <- eq
	}
	return(dat)
}
#'
# Run the function
dat <- behavior(envComb1 = c(-.05, -.05), envComb2= c(0.05, 0.05))
#+ fig.width = 5.5, fig.height = 4
# plot
par(family = 'serif', cex = 0.8, mar = c(4, 4, 1, 1))
plot(0,
     xlim = c(0, dim(dat)[1]),
     ylim = c(0, max(dat[,c(3:6)]) + 0.1),
     xlab = "time",
     ylab = "state proportion"
     )
lines(dat$B, col = 2, lwd = 1.5)
lines(dat$T, col = 3, lwd = 1.5)
lines(dat$M, col = 4, lwd = 1.5)
lines(dat$R, col = 5, lwd = 1.5)
legend(34, 0.25, c("B", "T", "M", "R"), col = 2:5, bty = "n", lty = 1)


#'#################################
#'# Effect of forest management on recovery resilience {#fm}
#'#################################
#'
#' Using a numerical approach, I will set forest management into the model by changing
#' the parameters related to each management process.
#' For example, __plantation__ can enhance colonization rate ($\alpha$), pre-commercial
#' __thinning__ can enhance competitive exclusion ($\theta$), __enrichment planting__ can enhance
#' succession rate ($\beta$), and __cutting__ can enhance disturbance rate ($\epsilon$; Figure 2).
#' Dotted line is the base natural process that occur without intervention.
#'
#' <center><img src="img/fig2.pdf" alt="Figure 2" style="width: 900px;"/></center>
#'
#' By simply increasing the value of each parameter related to management, we can simulate
#' the inclusion of management in the model processes (see figure above).
#' The __recovery resilience__, or time rate in which a system returns to equilibrium after
#' a disturbance, is measured by the largest real part of the __eigenvalue__.
#' The eigenvalue is optained by the Jacobian matrix and a nice example can be found in
#' this [vignette](https://cran.r-project.org/web/packages/rootSolve/vignettes/rootSolve.pdf)
#' of the `rootSolve` package.
#'
#' Update (October 2017)
#' ## Integrating forest management and biological mechanisms
#'
# function to plot parameter variation
managXmechan <- function(par1, ENV1, managVar, legend = NULL) { # managVar = vector with same lenght of parameters

  library(RColorBrewer)
  if(is.null(legend)) legend = "FALSE"

  # defining plot title (local environment)
  if(ENV1 < -0.9) {
    title <- "Boreal"
  }else title <- "Mixed"

  # parameter increasing
  dat <- data.frame(matrix(NA, ncol = 50, nrow = 9))
  dat$X1 = get_pars(ENV1 = ENV1, ENV2 = 0, params, int = 3)
  for(j in 1:9) {
    sq <- seq(0, managVar[j], length.out = 49)
    for(i in 1: length(sq)) {
      dat[j, i+1] <- dat[j, 1] + (dat[j, 1] * sq[i])
    }
  }

  # plot
  Col <- brewer.pal(4,"Dark2")

  par(family = "serif", mar = c(2, 4, 1, 1))
  plot(0, pch = "", xlim = c(0, 50), xaxt = "n", ylim = c(0, 1.8), type = "l", xlab = "", ylab = "parameter value")
  for(i in 1: length(par1)) points(c(1:50), dat[par1[i], ], type = "l", col = Col[i])
  mtext(title, 3)
  mtext("Foret management increasing (%) ", 1, cex = 1.1)
  mtext(expression(symbol("\256")), side = 1, line = 0, at = 45)
  if(legend == "TRUE") {
    legend("topleft", legend = c("alphaT", "betaT", "thetaT", "epsB"), lty = 1, col = Col, bty = "n", cex = 0.8)
  }
}
#'
# Used management
#+ fig.width = 10, fig.height = 4.5
par(mfrow = c(1, 2))
managXmechan(par1 = c(2, 4, 6, 7), ENV1 = - 1.54, managVar = c(1, 8, 0.4, 4, 0.4, 0.3, 3, 3, 3), legend = TRUE)
managXmechan(par1 = c(2, 4, 6, 7), ENV1 = .24, managVar = c(1, 1, 0.4, 0.4, 0.4, 0.3, 3, 3, 3))

# Extreme management
#+ fig.width = 10, fig.height = 4.5
par(mfrow = c(1, 2))
managXmechan(par1 = c(2, 4, 6, 7), ENV1 = - 1.54, managVar = c(1, 120, 0.4, 59, 0.6, 0.37, 42, 3, 3))
managXmechan(par1 = c(2, 4, 6, 7), ENV1 = .24, managVar = c(1, 0.3, 0.4, 1.7, 0.6, 0.34, 175, 3, 3))


# function to plot parameter variation
managXmechan <- function(par1, ENV1, managVar, legend = NULL) { # managVar = vector with same lenght of parameters

  library(RColorBrewer)
  if(is.null(legend)) legend = "FALSE"

  # defining plot title (local environment)
  if(ENV1 < -0.9) {
    title <- "Boreal"
  }else title <- "Mixed"

  ff <- function(x, y) {
    (0.02*x) + y
  }

  # parameter increasing
  dat <- data.frame(matrix(NA, ncol = 51, nrow = 9))
  dat$X1 = get_pars(ENV1 = ENV1, ENV2 = 0, params, int = 3)
  x <- 1:50

  for(j in 1:9) {
    sq <- seq(0, managVar[j], length.out = 49)
    for(i in 1: length(sq)) {
      dat[j, -1] <- ff(x, dat[j, 1])
    }
  }

  # plot
  Col <- brewer.pal(4,"Dark2")

  par(family = "serif", mar = c(2, 4, 1, 1))
  plot(0, pch = "", xlim = c(0, 50), xaxt = "n", ylim = c(0, 1.8), type = "l", xlab = "", ylab = "parameter value")
  for(i in 1: length(par1)) points(c(1:51), dat[par1[i], ], type = "l", col = Col[i])
  mtext(title, 3)
  mtext("Foret management increasing (%) ", 1, cex = 1.1)
  mtext(expression(symbol("\256")), side = 1, line = 0, at = 45)
  if(legend == "TRUE") {
    legend("topleft", legend = c("alphaT", "betaT", "thetaT", "epsB"), lty = 1, col = Col, bty = "n", cex = 0.8)
  }
}

#' The choice of parameters is delicate in this kind of "_sensitivity analysis_", then I chose
#' to do 2 different tests in the parametric variation to simulate the response of resilience
#' to the increasing in forest management.
#'
#' [`Test1`](#test1) varies a specific parameter (the tested one) from 0 to 1, and the
#' other parameters remains with the [original values](#ov) (fitted with field data).
#'
#' [`Test2`](#test2) varies a specific parameter from 0 to 1, and the other parameters
#' remains with a constant value `constPar`.
#' So we tested with four constant values: `0.2, 0.4, 0.6, 0.8`.
#'
#' ##Test 1 - paramenters constant with original value {#test1}

# Constant value for the main parameters from 0 to 1 + original value for the other parameters
# And environment 1 variation
int <- 3
parSeq <- seq(0.01, 1.1, 0.05)
Ev1 <- c(0, -1.35, -0.65, 0.45)
#'
# running eigenvalue to each parameter
eql <- as.list("NA")
df <- data.frame()
  for(k in 1: length(pars)) {
    for(i in 1: length(Ev1)) {
    	pars = get_pars(ENV1 = Ev1[i], ENV2 = 0, params, int = int)
    	for(j in 1: length(parSeq)) {
    		pars[k] = parSeq[j]
    		df[j, 1]	<- parSeq[j]
    		df[j, (i + 1)] <- get_eq(pars)$ev
    	}}
eql[[k]] <- df
}

#' **Simulations with ENV1 = 0 (Temperate for about Montreal latitude)**
#+ fig.width = 9.5, fig.height = 8
Pars <- c(expression(alpha), expression(alpha), expression(beta),
          expression(beta), expression(theta), expression(theta),
          expression(epsilon), expression(epsilon), expression(epsilon))
st <- c("B", "T", "B", "T", "", "T", "B", "T", "M")
par(family = 'serif', cex = 1.2, mfrow = c(3,3), mai = c(0.3, .5, .2, .2))
for(i in 1: length(pars)) {
	plot(eql[[i]]$V1, eql[[i]]$V2, type = "l", lwd = 1.5, xlab = "", ylab = "", ylim = c(-.45,0), xlim = c(0, 1.1))
  abline(v = get_pars(ENV1 = 0, ENV2 = 0, params, int = int), col = "gray55", lty = 3)
  if(i == 4) mtext(side = 2, "largest real part", line = 2.1, cex = 1.2)
  if(i == 2) mtext(Pars, at = get_pars(ENV1 = 0, ENV2 = 0, params, int = int), line = 0, cex = 0.75)
  legend("bottomleft", Pars[i], bty = "n")
  legend(0.005, -.41, st[i], bty = "n")
}
#'
#' **Simulations for different environments (-1.35 for Boreal, -0.65 for Mixed and 0.45 for Temperate)**
#+ fig.width = 9.5, fig.height = 8
par(family = 'serif', cex = 1.2, mfrow = c(3,3), mai = c(0.3, .5, .2, .2))
for(i in 1: length(pars)) {
	plot(eql[[i]]$V1, eql[[i]]$V3, type = "l", lwd = 1.5, xlab = "", ylab = "", ylim = c(-.5,0), xlim = c(0, 1.1))
  lines(eql[[i]]$V1, eql[[i]]$V4, lwd = 1.5, lty = 2)
  lines(eql[[i]]$V1, eql[[i]]$V5, lwd = 1.5, lty = 3)
  if(i == 4) mtext(side = 2, "largest real part", line = 2.1, cex = 1.2)
  if(i == 6) legend("bottomright", c("B (-1.35)", "M (-0.65)", "T (0.45)"), lty = 1:3, lwd = 1.5, bty = "n", cex = 1.2)
  legend("bottomleft", Pars[i], bty = "n")
  legend(0.01, -.45, st[i], bty = "n")
}
#'
#' ### Interaction between parameters

#' Some practices can change more than a parameter, so here I want to see the interaction
#' between these parameters influenced by the same practices.
#' The follow function can generate a plot to visualize the interaction between $\alpha$ Boreal
#' and $\alpha$ Temperate, $\beta$ Boreal and $\beta$ Temperate and $\theta$ and $\theta$ Temperate.
#'
# function to generate the plot interaction between parameters
parInt <- function(Par1, Par2, ENV1) {

    # parameter variation
    int <- 3
    par1 <- seq(0.01, 0.99, 0.05)
    par2 <- seq(0.01, 0.99, 0.05)

    # all possible combination
    BT <- expand.grid(p1 = par1, p2 = par2)

    # get eigenvalue for each combination
    for(i in 1:dim(BT)[1]) {
      pars = get_pars(ENV1 = ENV1, ENV2 = 0, params, int = int)
      pars[c(Par1, Par2)] <- c(BT$p1[i], BT$p2[i])
      BT[i, 3] <- get_eq(pars)$ev
    }

    # setting eigenvalue dependence color
    rbPal <- colorRampPalette(c('red','blue'))
    BT$Col <- rbPal(20)[as.numeric(cut(BT$V3,breaks = 15))]

    # plot
    library(grDevices)
    par(family = "serif", cex = 0.8, mar = c(4, 4, 1, 1))
    layout(matrix(1: 2, ncol = 2), width = c(2, 1), height = c(1, 1))
    plot(BT$p1, BT$p2, pch = 15, cex = 2.9, col = BT$Col, xlab = names(pars[Par1]), ylab = names(pars[Par2]))
    mtext(paste("ENV1 = ", ENV1, sep = ""), 3, cex = 1)

    # legend
    legend_image <- as.raster(matrix(rbPal(12), ncol=1))
    plot(c(0, 2),c(0, min(BT$V3)), type = 'n', axes = F, xlab = '', ylab = '', main = 'eigenvalue')
    text(x = 1, y = seq(max(BT$V3), min(BT$V3), l = 5), labels=round(seq(min(BT$V3), max(BT$V3),l=5),3))
    rasterImage(legend_image, min(BT$V3), min(BT$V3), 0, 0)
}
#'
#'**$\alpha$B and $\alpha$T interaction**
#'
#+ fig.width = 5, fig.height = 3
parInt(1, 2, ENV1 = -1.35)
parInt(1, 2, ENV1 = -.75)
parInt(1, 2, ENV1 = .5)
#'
#' **$\beta$B and $\beta$T interaction**
#'
#+ fig.width = 6, fig.height = 4
parInt(3, 4, ENV1 = -.75)
#'
#' **$\theta$ and $\theta$T interaction**
#+ fig.width = 5, fig.height = 3
parInt(5, 6, ENV1 = -1.35)
parInt(5, 6, ENV1 = -.75)
parInt(5, 6, ENV1 = .5)
#'
#' **$\epsilon$M and $\epsilon$T interaction**
#+ fig.width = 5, fig.height = 3
parInt(9, 8, ENV1 = -1.35)
parInt(9, 8, ENV1 = -.75)
parInt(9, 8, ENV1 = .5)
#'
#' ##Test 2 - parameters constant with the same value {#test2}

# Constant value for all parameters from 0 to 1
int <- 3
parSeq <- seq(0.01, 0.99, 0.05)
constPar <- c(0.2, 0.4, 0.6, 0.8) #Constant value of all other parameters
#'
# running eigenvalue to each parameter
pars = get_pars(ENV1 = 0, ENV2 = 0, params, int = int)
eql1 <- as.list("NA")
df <- data.frame()
for(k in 1: length(pars)) {
	pars = get_pars(ENV1 = 0, ENV2 = 0, params, int = int)
  for(l in 1: length(constPar)) {
    pars[-k] <- constPar[l]
  	for(j in 1: length(parSeq)) {
  		pars[k] = parSeq[j]
  		df[j, 1]	<- parSeq[j]
  		df[j, l + 1] <- get_eq(pars)$ev
  	 }
  }
eql1[[k]] <- df
}
#'

#+ fig.width = 9.5, fig.height = 8
# plot
source("/Users/wvieira/GitHub/NativeFunctions/addTransColor.R") #get transparance function
Pars <- c(expression(alpha), expression(alpha), expression(beta),
          expression(beta), expression(theta), expression(theta),
          expression(epsilon), expression(epsilon), expression(epsilon))
st <- c("B", "T", "B", "T", "", "T", "B", "T", "M")
par(family = 'serif', cex = 1.2, mfrow = c(3,3), mai = c(0.3, .5, .2, .2))
for(i in 1: length(pars)) {
	plot(0, type = "n", lwd = 1.5, xlab = "", ylab = "", xlim = c(0, 1), ylim = c(-0.35, 0.01))
  if(i == 4) mtext(side = 2, "largest real part", line = 2.1, cex = 1.2)
  if(i == 9) legend("bottomright", legend = constPar, lty = 1, col = c(1:4), bty = "n", cex = 1.2)
  legend("bottomleft", Pars[i], bty = "n")
  legend(0.01, -.32, st[i], bty = "n")
      points(eql1[[i]]$V1, eql1[[i]]$V2, type = "l", lwd = 1.5)
      points(eql1[[i]]$V1, eql1[[i]]$V3, type = "l", lwd = 1.5, col = 2)
      points(eql1[[i]]$V1, eql1[[i]]$V4, type = "l", lwd = 1.5, col = 3)
      points(eql1[[i]]$V1, eql1[[i]]$V5, type = "l", lwd = 1.5, col = 4)
    #add vertical lines
    for(j in 1:4) abline(v = constPar[j], col = addTrans(j, 120), lty = 3)
}

#'#################################
#'## Testing all parameters combination
#'#################################
#'
#' Now we saw the variation of each parameter, we want to see every parameter changing with
#' each other to test all combinations of parameters.
#' The sensitivity analysis helps us to test some of the following questions:
#'
#' - How is the relation between input and outputs of the model?
#'
#' - For what parameter are the outputs varying more?
#'
#' - Can we reduce uncertainty and hence increase robustness?
#'
#' All these questions will then be guided to measure the impact of forest management on
#' the parameters (based in [figure 2](#fm)).
#'
#' Using the `expand.grid` function, the parameters variation goes from 0.09 to 0.99:
#'

# parameters
parm <- seq(0.09, 0.99, 0.15)

# all possible combination
sim <- expand.grid(alphab = parm,
                   alphat = parm,
                   betab = parm,
                   betat = parm,
                   theta = parm,
                   thetat = parm,
                   eps = parm)
#'
dim(sim)

# get eigenvalue for each variation (time ~ 29 hrs in one core):
# system.time(for(i in 1: dim(sim)[1]) {
#   eq <- get_eq(sim[i, ])
#   sim[i, 8] <- eq$ev
#   sim[i, 9] <- eq$eq[1] #eq B
#   sim[i, 10] <- eq$eq[2] #eq T
#   sim[i, 11] <- eq$eq[3] #eq M
#   cat(100*i/dim(sim)[1], "%", "\r")
# })
#'
# save the output
# save(sim, file = "data/ev.RData")
#'
# load the output
load("data/evq.RData")

#'#################################
#'### Exploratory analysis
#'#################################
#'
#' I'm going to start with an analysis of variance type III to get what parameter
#' influenced more the eigenvalue.

library(car)

# Call anova eigenvalue (ev ~ all parameters variation + interactions)
aov <- Anova(lm1 <- lm(V8 ~
                        alphab + alphat + betab + betat + theta + thetat + eps +
                        alphab * alphat + alphab * eps + alphat * eps + betab * theta +
                        betab * thetat + betat * theta + betat * thetat,
                data = sim),
                singular.ok = FALSE,
                type = "III"
)

#'
#' This is
# table
knitr::kable(aov)
#'
#' This is
# function to calculate Omega²
source("/Users/wvieira/GitHub/NativeFunctions/CalcOmegaSq.R")
# get Omega²
om <- calcOmegaSq(aov, Order = TRUE)
#'
# table
knitr::kable(om)
#'
#'#################################
#'# Different management scenarios
#'#################################
#'
#' Here I will test different scenarios of management practices to understand their impact
#' on forest resilience.
#' Each scenario will be created following the same numerical approach explained in [this section](#fm).
#' Different from [test 1](#test1), where we consider the increase in management practices
#' for each parameter separated, here we will consider the interaction between practices
#' because a management practice is never alone (e.g. after harvesting an area, plantation
#' is expected; when thinning an area, future harvest and hence plantation is expected).
#'
#' Considering the four main management practices (plantation, thinning, Enrichment plantation,
#' and harvest), we will start with six different scenarios: (i) increasing harvest rate leads
#' to increase in plantation rate; (ii) increasing thinning rate, leads to increase harvest and
#' hence plantation rate. (iii) Increasing in enrichment leads to increase harvest and hence plantation,
#' (iv) increasing in thinning and enrichment rate leads to increase harvest and hence plantation,
#' (v) increasing thinning leads to increase harvest and (vi) just exploitation (harvest).
#'
#' Scenarios:
#'
#' 1. plantation + harvest
#' 2. plantation + thinning + harvest
#' 3. plantation + enrichment + harvest
#' 4. plantation + thinning + enrichment + harvest
#' 5. thinning + harvest
#' 6. harvest
#'
#' Each of these two scenarios, will be multiplied by two for management focusing in the (i)
#' Boreal forest or in the (ii) Temperate one.
#'
#' These scenarios are summarized in the following table:
#'
#' | Parameters | Scenario 1 | Scenario 1 | Scenario 2 | Scenario 2 | Scenario 3 | Scenario 3 | Scenario 4 | Scenario 4 | Scenario 5 | Scenario 5 | Scenario 6 | Scenario 6 | Scenario 6 |
#' |:----------:|:----------:|:----------:|:----------:|:----------:|:----------:|:----------:|:----------:|:----------:|:----------:|:----------:|:----------:|:----------:|:----------:|
#' |            |   **B**    |   **T**    |    **B**   |    **T**   |    **B**   |    **T**   |    **B**   |    **T**   |   **B**    |   **T**    |   **B**    |   **M**    |   **T**    |
#' |  $\alpha$B |  &#8593;   |     -      |  &#8593;   |     -      |  &#8593;   |     -      |  &#8593;   |     -      |            |            |            |            |            |
#' |  $\alpha$T |     -      |  &#8593;   |     -      |  &#8593;   |     -      |  &#8593;   |     -      |  &#8593;   |            |            |            |            |            |
#' |  $\beta$B  |     -      |            |            |            |  &#8593;   |     -      |  &#8593;   |     -      |            |            |            |            |            |
#' |  $\beta$T  |            |            |            |            |     -      |  &#8593;   |     -      |  &#8593;   |            |            |            |            |            |
#' |  $\theta$  |            |            |  &#8593;   |     -      |            |            |  &#8593;   |     -      |  &#8593;   |     -      |            |            |            |
#' |  $\theta$T |            |            |     -      |  &#8593;   |            |            |     -      |  &#8593;   |     -      |  &#8593;   |            |            |            |
#' | $\epsilon$B|  &#8593;   |     -      |  &#8593;   |     -      |  &#8593;   |     -      |  &#8593;   |     -      |  &#8593;   |     -      |  &#8593;   |     -      |     -      |
#' | $\epsilon$M|     -      |     -      |     -      |     -      |     -      |     -      |     -      |     -      |     -      |     -      |     -      |  &#8593;   |     -      |
#' | $\epsilon$T|     -      |  &#8593;   |     -      |  &#8593;   |     -      |  &#8593;   |     -      |  &#8593;   |     -      |  &#8593;   |     -      |     -      |  &#8593;   |
#'
#'
#' Using scenario 1 as an example, to increase harvest and plantation rate we will increase the parameters $\epsilon$
#' and $\alpha$, respectively. If the objective is to focus on Borel, we will increase $\alpha$B
#' and keep $\alpha$T with the original (or natural) value.
#'
#' **Model outputs**:
#'
#' Other than evaluationg the largest real part of the eigenvalue as a resilience measure,
#' I will evaluate two other outputs. First, the proportion of each state as a forest
#' productivity meansure. And second, the transition probability (from the transition matrix)
#' from T -> M and from M -> B. Using the transition probability I can estimate the migration
#' rate towards North, where a bigger transition means bigger adaptation to climate change.
#'
#' ## Relative variation: parameter variation using the original values
#'
#' Varying the parameters relatively is a good way to take the original value of each parameter
#' into account. However, this variation may not represent very well the effect of
#' management practice (e.g. plantation can have a bigger impact than harvest because $\alpha$
#' is much bigger `alpha = 0.9` than $\epsilon$ `eps = 0.09`. For that, I'm going to increase the impact
#' of harvest on disturbance because the value of $\epsilon$ is too small comparing with other
#' parameters.
#'
# Sequence from 0 to 100% of Increase
sq <- seq(0.1, 1, 0.1) # 0 to 100% for all parameters except epsilon
sqH <- seq(0.3, 3, 0.3) # 30 to 300% for epsilon
#'
# Varing temperature (ENV1) from -2 to 2 and then plot each parameter variation
tempVar <- seq(-2, 2, length.out = 100)
dat.Par <- data.frame(matrix(NA, ncol = length(pars) - 2, nrow = length(tempVar)))
names(dat.Par) <- names(pars)[1:7]

#loop
for(i in 1: length(tempVar)) {
  dat.Par[i, ] = get_pars(ENV1 = tempVar[i], ENV2 = 0, params, int = 3)[1:7]
}

# plot
#+ fig.width = 6, fig.height = 5.5
plot(tempVar, dat.Par[, 1], type = "l", lwd = 1.5,
     ylab = "parameter value", xlab = "Temperature variation", ylim = c(0, 1))
for(i in 2: 7) {
  points(tempVar, dat.Par[, i], type = "l", col = i, lwd = 1.5)
}
legend("topright", names(pars)[1:7], col = 1:7, lwd = 1.5, bty = "n")
abline(v = c(-1.35, -0.65, 0.45), lty = 3, col = "gray")
mtext(c("Boreal", "Mixed", "Temperate"), 3, at = c(-1.35, -0.65, 0.45))

# define relative variation for each parameter
managSim <- function(Par1, ENV1, managVar) {

  # defining management increasing for each parameter
  oP <- 0
  managV <- rep(NA, length(pars))
  managV[Par1] <- managVar
  managV[-Par1] <- 0

  # defining plot title (local environment)
  if(ENV1 < -0.9) {
    title <- "Boreal"
  }else if(ENV1 > -0.2) {
    title <- "Temperate"
  }else title <- "Mixed"

  # data frame
  dat <- data.frame(matrix(NA, ncol = 16, nrow = length(pars)))
  dat$X1 = get_pars(ENV1 = ENV1, ENV2 = 0, params, int = 3)
  for(j in 1: length(pars)) {
    sq <- seq(0.1, managV[j], length.out = 15)
    for(i in 1: length(sq)) {
      dat[j, i+1] <- dat[j, 1] + (dat[j, 1] * sq[i])
    }
  }

  # data frame
  # X2 = Par1 eigenvalue; X3:X6 = Par1 state proportion
  egv <- data.frame(matrix(ncol = 6, nrow = dim(dat)[2]))
  egv$X1 <- seq(1, dim(dat)[2])
  pars <- get_pars(ENV1 = ENV1, ENV2 = 0, params, int = 3)

  # get eigenvalue and satate proportion for Par1
  for(i in 1: (length(sq) + 1)) {
    pars[Par1] <- dat[Par1, i]
    eq <- get_eq(pars)
    egv[i, 2] <- eq$ev #eigenvalue
    egv[i, c(3:5)] <- eq$eq # B, T, M
    egv[i, 6] <- 1 - sum(eq$eq) # R
  }

  # get probability transition
  prob <- data.frame(matrix(ncol = 3, nrow = length(sq) + 1))
  prob[1] <- seq(1, dim(dat)[2])
  pars <- get_pars(ENV1 = ENV1, ENV2 = 0, params, int = 3)
  eq <- get_eq(pars)[[1]]
  for(k in 1: (length(sq)+ 1)) {
    pars[Par1] <- dat[Par1, k]
    eq <- get_eq(pars)[[1]]
    MAT <- get_matrix(ENV1 = ENV1, ENV2 = 0, pars = pars, eq = eq)$MAT
    prob[k, 2] <- MAT[2, 3]
    prob[k, 3] <- MAT[3, 4]
  }

  # plot 1
  par(mfrow = c(1, 3), family = "serif", mar = c(2, 4, 1, 1), cex = 1.2)
  plot(egv$X1, egv$X2, type = "l", xaxt = "n", ylim = c(-0.27, -0.01), lwd = 1.75,
       xlab = "", ylab = "Largest real part")

  #plot 2
  plot(egv$X1, egv$X3, type = "l", ylim = c(0, 1), lwd = 1.75, xaxt = "n",
       xlab = "", ylab = "Occupancy")
  lines(egv$X1, egv$X4, lwd = 1.75, col = 2)
  lines(egv$X1, egv$X5, lwd = 1.75, col = 3)
  lines(egv$X1, egv$X6, lwd = 1.75, col = 4)

  legend("topleft", c("B", "T", "M", "R"), lty = 1, col = c(1:4), bty = "n", cex = 0.9)
  mtext(paste(title, " (", ENV1, ")", sep = ""), 3)
  mtext("Foret management increasing (%)", 1, cex = 1.1)

  #plot 3
  plot(prob$X1, prob$X2, xaxt = "n", type = "l", lwd = 1.75, ylim = c(0, 0.14),
       col = 3, ylab = "Transition probability", xlab = "")
  lines(prob$X1, prob$X3, lwd = 1.75, col = 2)
  legend("topleft", c("B -> M", "M -> T"), lty = 1, col = c(3, 2), bty = "n", cex = 0.9)
}
#'
#' ### Different scenarios
#'
#' #### Scenario 1: &#8593; plantation and &#8593; harvest
#'
#+ fig.width = 12, fig.height = 3.8
managSim(Par1 = c(2, 7), ENV1 = -1.54, managVar = c(1, 3))
managSim(Par1 = c(2, 7, 9), ENV1 = 0.24, managVar = c(1, 3, 3))
#'
#' #### Scenario 2: &#8593; plantation, &#8593; thinning (M) and &#8593; harvest
#+ fig.width = 12, fig.height = 3.8
managSim(Par1 = c(2, 6, 7), ENV1 = -1.54, managVar = c(1, 0.5, 3))
managSim(Par1 = c(2, 6, 9), ENV1 = 0.24, managVar = c(1, 0.5, 3))
#'
#' #### Scenario 3: &#8593; plantation; &#8593; enrichment (T) and &#8593; harvest
#+ fig.width = 12, fig.height = 3.8
managSim(Par1 = c(2, 4, 7), ENV1 = -1.54, managVar = c(1, 1, 3))
managSim(Par1 = c(2, 4, 9), ENV1 = 0.24, managVar = c(1, 1, 3))
#'
#' #### Scenario 4: &#8593; plantation; &#8593; thinning (M); &#8593; enrichment (T) and &#8593; harvest
#+ fig.width = 12, fig.height = 3.8
managSim(Par1 = c(2, 6, 4, 7), ENV1 = -1.54, managVar = c(1, 0.6, 0.6, 3))
managSim(Par1 = c(2, 6, 4, 9), ENV1 = 0.24, managVar = c(1, 0.6, 0.6, 3))
#'
#' #### Scenario 5: &#8593; thinning (M); harvest
#+ fig.width = 12, fig.height = 3.8
managSim(Par1 = c(6, 7), ENV1 = -1.54, managVar = c(0.6, 3))
managSim(Par1 = c(6, 9), ENV1 = 0.24, managVar = c(0.6, 3))
#'
#' #### Scenario 6: &#8593; harvest
#+ fig.width = 12, fig.height = 3.8
managSim(Par1 = 7, ENV1 = -1.54, managVar = 3)
managSim(Par1 = c(7, 9), ENV1 = -1, managVar = c(0 ,0))
#'
#' ## Impact of each practice and each management scenario on resilience
#' '(figure to be used in the manuscript)
#'
# setting colors (colorFriendly = TRUE)
library(RColorBrewer)
colPrac <- brewer.pal(4,"Dark2")
colManag <- brewer.pal(6, "Set1")

# function
resilManag <- function(scenario, ENV1, managVar, legend = NULL) { # scenario = list; managVar = vector with same lenght of parameters

  # defining plot title (local environment)
  if(ENV1 < -0.9) {
    title <- "Boreal"
  }else title <- "Mixed"

  if(is.null(legend)) legend <- "FALSE"

  # parameter increasing
  dat <- data.frame(matrix(NA, ncol = 50, nrow = 9))
  dat$X1 = get_pars(ENV1 = ENV1, ENV2 = 0, params, int = 3)
  for(j in 1:9) {
    sq <- seq(0.1, managVar[j], length.out = 49)
    for(i in 1: length(sq)) {
      dat[j, i+1] <- dat[j, 1] + (dat[j, 1] * sq[i])
    }
  }

  # table for output
  egv <- data.frame(matrix(ncol = length(scenario) + 1, nrow = dim(dat)[2]))
  egv$X1 <- seq(1, dim(dat)[2])

  for(w in 1: length(scenario)) {

    # define parameters for each scenario
    Par1 <- scenario[[w]]

    # model running
    pars <- get_pars(ENV1 = ENV1, ENV2 = 0, params, int = 3)
    # get eigenvalue and satate proportion for Par1
    for(i in 1: dim(dat)[2]) {
      pars[Par1] <- dat[Par1, i]
      eq <- get_eq(pars)
      egv[i, w + 1] <- eq$ev #eigenvalue
    }
  }

  # define color
  if(length(scenario) == 4) {
    color <- colPrac
  }else color <- colManag

  # plot
  par(family = "serif", mar = c(2, 4, 1, 1), cex = 1.2)
  plot(egv$X1, egv$X2, type = "l", col = color[1], xaxt = "n", ylim = c(-.135, 0), lwd = 1.75,
       xlab = "", ylab = "Largest real part")
  for(k in 3: (length(scenario) + 1)) points(egv$X1, egv[,k], type = "l", lty = k-1, lwd = 1.75, col = color[k - 1])
  if(legend == TRUE) {
    if(length(scenario) == 4) {
    legend("bottomleft", c("Plantation", "Enrichment", "Thinning", "Harvest"), lty = 1:4, col = color, bty = "n", cex = 0.75)
  }else legend("bottomleft", c("P + H", "P + T + H", "P + E + H", "P + T + E + H", "T + H", "H"), lty = 1:6, col = color, bty = "n", cex = 0.75)
  }
  mtext(title, 3)
  mtext("Foret management increasing (%) ", 1, cex = 1.1)
  mtext(expression(symbol("\256")), side = 1, line = 0, at = 45)

}
#'

managVar <- c(1, 2, .4, .4, .55, .55, 3, 3, 3)

# plot 1
scenario1 = as.list(c(2, 4, 6, 7))
scenario2 = as.list(c(2, 4, 6, 9))
#+ fig.width = 10, fig.height = 4.5
par(mfrow = c(1, 2))
resilManag(scenario1, ENV1 = -1.54, managVar, legend = TRUE)
mtext(expression(symbol("\255")), side = 3, line = -2.5, at = 27.5, cex = 1.8)
resilManag(scenario2, ENV1 = 0.24, managVar)
mtext(expression(symbol("\255")), side = 3, line = -2.5, at = 11, cex = 1.8)

# plot 2
scenarioB <- list(c(2, 7), c(2, 6, 7), c(2, 4, 7), c(2, 4, 6, 7), c(6, 7), 7)
scenarioM <- list(c(2, 9), c(2, 6, 9), c(2, 4, 9), c(2, 4, 6, 9), c(6, 9), 9)

#+ fig.width = 10, fig.height = 4.5
par(mfrow = c(1, 2))
resilManag(scenarioB, ENV1 = -1.54, managVar, legend = TRUE)
mtext(expression(symbol("\255")), side = 3, line = -2.5, at = 24, cex = 1.8)
mtext(expression(symbol("\255")), side = 3, line = -2.5, at = 28, cex = 1.8)
resilManag(scenarioM, ENV1 = 0.24, managVar)
mtext(expression(symbol("\255")), side = 3, line = -2.5, at = 9, cex = 1.8)

#######################################################################

# function
resilManag <- function(scenario, ENV1, managVar1 = NULL, managVar2 = NULL, incrMethod, legend = NULL) { # scenario = list; managVar = vector with same lenght of parameters

  # setting colors (colorFriendly = TRUE)
  library(RColorBrewer)
  colPrac <- brewer.pal(4,"Dark2")
  colManag <- brewer.pal(6, "Set1")

  # defining plot title (local environment)
  if(ENV1 < - 1.2) {
    title <- "Boreal"
  }else if(ENV1 < 0.4) {
    title <- "Mixed"
  }else title <- "Temperate"

  if(is.null(legend)) legend <- "FALSE"

  # parameter increasing (two methods, first increasing by % and second linearly)
    # method 1
  if(incrMethod == 1) {

    ff <- function(x, y) {
      (managVar1*x) + y
    }

    dat <- data.frame(matrix(NA, ncol = 51, nrow = 9))
    dat$X1 = get_pars(ENV1 = ENV1, ENV2 = 0, params, int = 3)
    x <- 1:50

    for(j in 1:9) {
      for(i in 1: 49) {
        dat[j, -1] <- ff(x, dat[j, 1])
      }
    }
  }
    # method 2
  if(incrMethod == 2) {
    dat <- data.frame(matrix(NA, ncol = 51, nrow = 9))
    dat$X1 = get_pars(ENV1 = ENV1, ENV2 = 0, params, int = 3)
    for(j in 1:9) {
      sq <- seq(0.1, managVar2[j], length.out = 50)
      for(i in 1: length(sq)) {
        dat[j, i+1] <- dat[j, 1] + (dat[j, 1] * sq[i])
      }
    }
  }

  # table for output
  egv <- data.frame(matrix(ncol = length(scenario) + 1, nrow = dim(dat)[2]))
  egv$X1 <- seq(1, dim(dat)[2])

  for(w in 1: length(scenario)) {

    # define parameters for each scenario
    Par1 <- scenario[[w]]

    # model running
    pars <- get_pars(ENV1 = ENV1, ENV2 = 0, params, int = 3)
    # get eigenvalue and satate proportion for Par1
    for(i in 1: dim(dat)[2]) {
      pars[Par1] <- dat[Par1, i]
      eq <- get_eq(pars)
      egv[i, w + 1] <- eq$ev #eigenvalue
    }
  }

  # define color
  if(length(scenario) == 4) {
    color <- colPrac
  }else color <- colManag

  # plot
  par(mar = c(2, 3, 0.8, 0.5), mgp = c(1.5, 0.3, 0), tck = -.01, family = 'serif', cex = 1.2)
  plot(egv$X1, egv$X2, type = "l", col = color[1], xaxt = "n", ylim = c(-.135, 0), lwd = 1.75,
       xlab = "", ylab = "Largest real part")
  for(k in 3: (length(scenario) + 1)) points(egv$X1, egv[,k], type = "l", lty = k-1, lwd = 1.75, col = color[k - 1])
  if(legend == TRUE) {
    if(length(scenario) == 4) {
    legend("bottomleft", c("Plantation", "Enrichment", "Thinning", "Harvest"), lty = 1:4, col = color, bty = "n", cex = 0.8)
    }else legend("bottomleft", c("P + H", "P + T + H", "P + E + H", "P + T + E + H", "T + H", "H"), lty = 1:6, col = color, bty = "n", cex = 0.8)
  }
  mtext(title, 3)
  mtext("Foret management increasing (%) ", 1, cex = 1.1)
  mtext(expression(symbol("\256")), side = 1, line = 0, at = 45)

}
#'
#' With the last function, we have two ways to calculate the increase in the parameters to simulate
#' the forest management. The first way is a linear increase with a constant slop determined by `managVar1`.
#' The second way is a proportional increase with the rate of increase varing from 50% (for thinning)
#' to 300% (for harvest). The following plots show these parameter variation.

# method 1
ff <- function(x, y) {
  (managVar1*x) + y
}

managVar1 = 0.005
dat <- data.frame(matrix(NA, ncol = 51, nrow = 9))
dat$X1 = get_pars(ENV1 = -1, ENV2 = 0, params, int = 3)
x <- 1:50

for(j in 1:9) {
  for(i in 1: 49) {
    dat[j, -1] <- ff(x, dat[j, 1])
  }
}

# method 2
managVar2 <- c(1, 1, .55, .55, .55, .55, 3, 3, 3)
dat1 <- data.frame(matrix(NA, ncol = 51, nrow = 9))
dat1$X1 = get_pars(ENV1 = -1, ENV2 = 0, params, int = 3)
for(j in 1:9) {
  sq <- seq(0.1, managVar2[j], length.out = 50)
  for(i in 1: length(sq)) {
    dat1[j, i+1] <- dat1[j, 1] + (dat1[j, 1] * sq[i])
  }
}
#'
# plot
#+ fig.width = 10, fig.height = 4.5
par(mfrow = c(1, 2), mar = c(2.8, 3, 1, 0.5), mgp = c(1.5, 0.3, 0), tck = -.01, family = 'serif', cex = 1.2)
plot(1:51, dat[1, ], type = "l", ylim = c(min(dat), max(dat)), xlab = "Forest management")
for(i in 2:7) points(1:51, dat[i, ], type = "l", col = i)
mtext("Method 1", 3)
plot(1:51, dat1[1, ], type = "l", ylim = c(min(dat1), max(dat1)), xlab = "Forest management",)
for(i in 2:7) points(1:51, dat1[i, ], type = "l", col = i)
mtext("Method 2", 3)
legend("topleft", names(pars)[1:7], col = 1:7, lwd = 1.5, bty = "n", cex = 0.9)
#'
# Pratices
  # method 1
#+ fig.width = 10, fig.height = 9
  par(mfrow = c(2, 2))
  managVar1 = 0.005
  resilManag(as.list(c(2, 4, 6, 7)), ENV1 = -1.54, managVar1, incrMethod = 1, legend = TRUE)
  resilManag(as.list(c(2, 4, 6, 7)), ENV1 = -1.09, managVar1, incrMethod = 1)
  resilManag(as.list(c(2, 4, 6, 9)), ENV1 = 0.24, managVar1, incrMethod = 1)
  resilManag(as.list(c(2, 4, 6, 9)), ENV1 = 0.53, managVar1, incrMethod = 1)
  # method 2
  par(mfrow = c(2, 2))
  managVar2 <- c(1, 1, .55, .55, .55, .55, 3, 3, 3)
  resilManag(as.list(c(2, 4, 6, 7)), ENV1 = -1.54, managVar2 = managVar2, incrMethod = 2, legend = TRUE)
  resilManag(as.list(c(2, 4, 6, 7)), ENV1 = -1.09, managVar2 = managVar2, incrMethod = 2)
  resilManag(as.list(c(2, 4, 6, 9)), ENV1 = 0.24, managVar2 = managVar2, incrMethod = 2)
  resilManag(as.list(c(2, 4, 6, 9)), ENV1 = 0.53, managVar2 = managVar2, incrMethod = 2)

# Management interaction
  # method 1
  scenarioB <- list(c(2, 7), c(2, 6, 7), c(2, 4, 7), c(2, 4, 6, 7), c(6, 7), 7)
  scenarioM <- list(c(2, 9), c(2, 6, 9), c(2, 4, 9), c(2, 4, 6, 9), c(6, 9), 9)
  par(mfrow = c(2, 2))
  managVar1 = 0.005
  resilManag(scenarioB, ENV1 = -1.54, managVar1, incrMethod = 1, legend = TRUE)
  resilManag(scenarioB, ENV1 = -1.09, managVar1, incrMethod = 1)
  resilManag(scenarioM, ENV1 = 0.24, managVar1, incrMethod = 1)
  resilManag(scenarioM, ENV1 = 0.53, managVar1, incrMethod = 1)
  # method 2
  par(mfrow = c(2, 2))
  managVar2 <- c(1, 1, .55, .55, .55, .55, 3, 3, 3)
  resilManag(scenarioB, ENV1 = -1.54, managVar2 = managVar2, incrMethod = 2, legend = TRUE)
  resilManag(scenarioB, ENV1 = -1.09, managVar2 = managVar2, incrMethod = 2)
  resilManag(scenarioM, ENV1 = 0.24, managVar2 = managVar2, incrMethod = 2)
  resilManag(scenarioM, ENV1 = 0.53, managVar2 = managVar2, incrMethod = 2)
