#'---
#'title: "Linking SDMs and forest management"
#'author: Willian Vieira
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
#' <center><img src="img/fig1.png" alt="Figure 1" style="width: 500px;"/></center>
#'
#' Each parameter represents a different process: $\alpha$ is the colonisation
#' rate (R -> B, T, M), $\beta$ the succession (B, T -> M), $\theta$ the exclusion
#' extintion (M -> B, T) and $\epsilon$ the pertubation (B, T, M -> R).
#'
#'# Example of running model
#'
# call the model
source("vissault_model.R")
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
transition_matrix = get_matrix(ENV1 = 0, ENV2 = 0, params, int = 3)
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

# loop
for(i in 2:length(time)) {
  eq = runsteady(y = y, func = model, parms = pars, times = c(0, i))[[1]]
  eq[4] <- 1 - eq[1] - eq[2] - eq[3]
  dm2[i, ] <- eq
}

#+ fig.width = 6.2, fig.height = 4.5
# plot
par(family = 'serif', cex = 0.8)
plot(time, dm2[, 1],
     type = "l",
     ylim = c(0, max(dm2)),
     col = 2,
     ylab = "State proportion")
lines(dm2[, 2], col = 3)
lines(dm2[, 3], col = 4)
lines(dm2[, 4], col = 5)
legend(44, 0.18, c("B", "T", "M", "R"), col = 2:5, lty = 1, bty = "n", cex = 0.8)

#'
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
#+ fig.width = 6.2, fig.height = 4.5
# plot
par(family = 'serif', cex = 0.8)
plot(0,
     xlim = c(0, dim(dat)[1]),
     ylim = c(0, max(dat[,c(3:6)]) + 0.1),
     xlab = "time",
     ylab = "state proportion"
     )
lines(dat$B, col = 2)
lines(dat$T, col = 3)
lines(dat$M, col = 4)
lines(dat$R, col = 5)
legend(34, 0.25, c("B", "T", "M", "R"), col = 2:5, bty = "n", lty = 1)


#'#################################
#'# Effect of forest management on resilience {#fm}
#'#################################
#'
#' Using a numerical approach, I will set forest management into the model by changing
#' the parameters related to each management process.
#' For example, __plantation__ can enhance colonization ($\alpha$), pre-commercial
#' __thinning__ can enhance both competitive exclusion ($\theta$) and succession ($\beta$),
#' and __cutting__ can enhance disturbance ($\epsilon$; Figure 2).
#' Dotted line is the base natural process that occur without intervention.
#'
#' <center><img src="img/fig2.pdf" alt="Figure 2" style="width: 700px;"/></center>
#'
#' By simply increasing the value of each parameter related to management, we can simulate
#' the inclusion of management in the model processes (see figure above).
#' The __recovery resilience__, or time rate in which a system returns to equilibrium after
#' a disturbance, is measured by the largest real part of the __eigenvalue__.
#' The eigenvalue is optained by the Jacobian matrix and a nice example can be found in
#' this [vignette](https://cran.r-project.org/web/packages/rootSolve/vignettes/rootSolve.pdf)
#' of the `rootSolve` package.
#'
#' The choice of parameters is delicate in this kind of "_sensitivity analysis_", then I chose
#' to do 2 different tests in the parametric variation to simulate the response of resilience
#' to the increasing in forest management.
#'
#' [`Test1`](#test1) varies a specific parameter (the tested one) from 0 to 1, and the
#' other parameters remains with the [original values](#ov) (fitted with field data).
#'
#' [`Test2`](#test2) varies a specific parameter from 0 to 1, and the other parameters
#' remains with a constant value `constPar`.
#' So we tested with four constant values: `0.2, 0.5, 0.8, 1.2`.
#'
#' ###Test 1 {#test1}

# Constant value for the main parameters from 0 to 1 + original value for the other parameters
int <- 2
parSeq <- seq(0.01, 1.1, 0.05)
#'
# running eigenvalue to each parameter
pars = get_pars(ENV1 = 0, ENV2 = 0, params, int = int)
eql <- as.list("NA")
df <- data.frame()
for(k in 1: length(pars)) {
	pars = get_pars(ENV1 = 0, ENV2 = 0, params, int = int)
	for(j in 1: length(parSeq)) {
		pars[k] = parSeq[j]
		df[j, 1]	<- parSeq[j]
		df[j, 2] <- get_eq(pars)$ev
	}
eql[[k]] <- df
}

#+ fig.width = 9.5, fig.height = 8
# plot
Pars <- c(expression(alpha), expression(alpha), expression(beta),
          expression(beta), expression(theta), expression(theta),
          expression(epsilon))
st <- c("B", "T", "B", "T", "", "T", "")
par(family = 'serif', cex = 1.2, mfrow = c(3,3), mai = c(0.3, .5, .2, .2))
for(i in 1:7) {
	plot(eql[[i]], type = "l", lwd = 1.5, xlab = "", ylab = "", ylim = c(-.31,0), xlim = c(0, 1.1))
  abline(v = get_pars(ENV1 = 0, ENV2 = 0, params, int = int), col = "gray", lty = 3, cex = 0.8)
  if(i == 4) mtext(side = 2, "largest real part", line = 2.1, cex = 1.2)
  legend("bottomleft", Pars[i], bty = "n")
  legend(0.01, -.29, st[i], bty = "n")
}
plot(c(0, 1.1), c(-.31,0), ann = FALSE, axes = FALSE, type = "n")
abline(v = get_pars(ENV1 = 0, ENV2 = 0, params, int = int), col = "gray", lty = 3, cex = 0.8)
text(x = get_pars(ENV1 = 0, ENV2 = 0, params, int = int), y = -0.1, Pars)

#' ###Test 2 {#test2}

# Constant value for all parameters from 0 to 1
int <- 2
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
          expression(epsilon))
st <- c("B", "T", "B", "T", "", "T", "")
par(family = 'serif', cex = 1.2, mfrow = c(3,3), mai = c(0.3, .5, .2, .2))
for(i in 1:7) {
	plot(0, type = "n", lwd = 1.5, xlab = "", ylab = "", xlim = c(0, 1), ylim = c(-0.35, 0.01))
  if(i == 4) mtext(side = 2, "largest real part", line = 2.1, cex = 1.2)
  legend("bottomleft", Pars[i], bty = "n")
  legend(0.01, -.32, st[i], bty = "n")
      points(eql1[[i]]$V1, eql1[[i]]$V2, type = "l", lwd = 1.5)
      points(eql1[[i]]$V1, eql1[[i]]$V3, type = "l", lwd = 1.5, col = 2)
      points(eql1[[i]]$V1, eql1[[i]]$V4, type = "l", lwd = 1.5, col = 3)
      points(eql1[[i]]$V1, eql1[[i]]$V5, type = "l", lwd = 1.5, col = 4)
    #add vertical lines
    for(j in 1:4) abline(v = constPar[j], col = addTrans(j, 120), lty = 3)
}
#add second legend
plot(c(-1, 1), c(-1, 1), ann = FALSE, axes = FALSE, type = "n")
legend(-1, 1, constPar, lty = 1, col = c(1:4), bty = "n", cex = 1.2)

#'
#'## Testing all parameters combination
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

# expand grid function
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

#' ## Exploratory analysis
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

# Call anova state B + M + T (ev ~ all parameters variation + interactions)
aov1 <- Anova(lm1 <- lm(V9 + V10 + V11 ~
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
knitr::kable(cbind(aov[c(2,3)],aov1[3]))
#'
#' This is
# function to calculate Omega²
source("/Users/wvieira/GitHub/NativeFunctions/CalcOmegaSq.R")
# get Omega²
om <- calcOmegaSq(aov)
om1 <- calcOmegaSq(aov1)
#'
# table
knitr::kable(cbind(om[c(1,2)],om1[2]))
