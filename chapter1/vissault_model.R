#################################
# Define the model
#################################
model = function(t, y, params) {
	with(as.list(c(y, params)), {

		# Fraction of empty patches converted into the different states following a disturbance
		pB = alphab * (B + M)
		pT = alphat * (T + M)
		pM = pB * pT
		pB_ = pB * (1 - pT) #what is it exaclty?
		pT_ = pT * (1 - pB) #what is it exaclty?

		# Regeneration state
		R = 1 - B - T - M

		# Differential equations describing the dynamics of the state variables
		dBdt = pB_ * R + theta * (1 - thetat) * M - betat * (T + M) * B - eps * B
		dTdt = pT_ * R + theta * thetat * M - betab * (B + M) * T - eps * T
		dMdt = pM * R + betab * (B + M) * T + betat * (T + M) * B - theta * M - eps * M
		list(c(dBdt, dTdt, dMdt))
		})
	}

#################################
# Local stability
#################################
get_eq = function(params, y = NULL) {

	library(rootSolve)

	# Initial conditions
	if(is.null(y)) {
		y = c(B = 0.25, T = 0.25, M = 0.25)
	}else(y = y)

	# Get the equilibrium
	eq = runsteady(y = y, func = model, parms = params, times = c(0, 1000))[[1]]

	# Compute the Jacobian
	J = jacobian.full(y = eq, func = model, parms = params)

	# Stability
	ev = max(Re(eigen(J)$values)) #in case of complex eigenvalue, using Re to get the first real part

	return(list(eq = eq, ev = ev))
}

#################################
# Wrapper to collect parameters for a given set of environmental conditions
#################################
logit_reverse <- function(x) exp(x) / (1 + exp(x))

get_pars = function(ENV1, ENV2, params, int) {

	logit_alphab 	= params["ab0", 1] + params["ab1", 1] * ENV1 + params["ab2", 1] * ENV2 + params["ab3", 1] * ENV1^2 + params["ab4",1]*ENV2^2 + params["ab5",1]*ENV1^3 + params["ab6",1]*ENV2^3
	logit_alphat 	= params["at0", 1] + params["at1", 1] * ENV1 + params["at2", 1] * ENV2 + params["at3", 1] * ENV1^2 + params["at4",1]*ENV2^2 + params["at5",1]*ENV1^3 + params["at6",1]*ENV2^3
	logit_betab 	= params["bb0", 1] + params["bb1", 1] * ENV1 + params["bb2", 1] * ENV2 + params["bb3", 1] * ENV1^2 + params["bb4",1]*ENV2^2 + params["bb5",1]*ENV1^3 + params["bb6",1]*ENV2^3
	logit_betat 	= params["bt0", 1] + params["bt1", 1] * ENV1 + params["bt2", 1] * ENV2 + params["bt3", 1] * ENV1^2 + params["bt4",1]*ENV2^2 + params["bt5",1]*ENV1^3 + params["bt6",1]*ENV2^3
	logit_theta		= params["th0", 1] + params["th1", 1] * ENV1 + params["th2", 1] * ENV2 + params["th3", 1] * ENV1^2 + params["th4",1]*ENV2^2 + params["th5",1]*ENV1^3 + params["th6",1]*ENV2^3
	logit_thetat	= params["tt0", 1] + params["tt1", 1] * ENV1 + params["tt2", 1] * ENV2 + params["tt3", 1] * ENV1^2 + params["tt4",1]*ENV2^2 + params["tt5",1]*ENV1^3 + params["tt6",1]*ENV2^3
	logit_eps 		= params["e0", 1]  + params["e1", 1]  * ENV1 + params["e2", 1]  * ENV2 + params["e3", 1]  * ENV1^2 + params["e4",1] *ENV2^2 + params["e5",1] *ENV1^3 + params["e6",1] *ENV2^3

	alphab = 1-(1-logit_reverse(logit_alphab))^int
	alphat = 1-(1-logit_reverse(logit_alphat))^int
	betab  = 1-(1-logit_reverse(logit_betab))^int
	betat  = 1-(1-logit_reverse(logit_betat))^int
	theta  = 1-(1-logit_reverse(logit_theta))^int
	thetat = 1-(1-logit_reverse(logit_thetat))^int
	eps    = 1-(1-logit_reverse(logit_eps))^int

	return(c(alphab = alphab, alphat = alphat, betab = betab, betat = betat, theta = theta, thetat = thetat, eps = eps))

}

#################################
# COLLECT THE TRANSITION MATRIX
#################################

	get_matrix = function(ENV1, ENV2, params, int) {

		pars = get_pars(ENV1,ENV2,params,int)
		eq = get_eq(pars)[[1]]

		MAT = matrix(nr = 4, nc = 4)

		B = eq[1]
		T = eq[2]
		M = eq[3]
		R = 1 - B - M - T
		names(R) <- "R" #Renaming R because it would get "B" name instanted

		pB = pars["alphab"]*(B+M)
		pT = pars["alphat"]*(T+M)
		pM = pB*pT
		pB_ = pB*(1-pT)
		pT_ = pT*(1-pB)

		# ORDER IN THE MATRIX: R, B, M, T

		MAT[1,1] = 1-(pB_*R + pM*R + pT_*R)
		MAT[1,2] = pB_*R
		MAT[1,3] = pM*R
		MAT[1,4] = pT_*R

		MAT[2,1] = pars["eps"]*B
		MAT[2,2] = 1 - pars["eps"]*B - pars["betab"]*(B+M)*T
		MAT[2,3] = pars["betab"]*(B+M)*T
		MAT[2,4] = 0

		MAT[3,1] = pars["eps"]*M
		MAT[3,2] = pars["theta"]*(1-pars["thetat"])*M
		MAT[3,3] = 1 - pars["eps"]*M - pars["theta"]*M
		MAT[3,4] = pars["theta"]*(pars["thetat"])*M

		MAT[4,1] = pars["eps"]*T
		MAT[4,2] = 0
		MAT[4,3] = pars["betat"]*(T+M)*B
		MAT[4,4] = 1 - pars["eps"]*T - pars["betat"]*(T+M)*B

		# Rename Matrix col and row
		nm <- c("R", "B", "M", "T")
		rownames(MAT) <- nm
		colnames(MAT) <- nm

		return(list(eq=c(R,B,M,T),MAT=MAT))

	}

#################################
# RUN THE MODEL TO EQUILIBRIUM
#################################

params = read.table("pars.txt", row.names = 1)
pars = get_pars(ENV1 = 0, ENV2 = 0, params, int = 3)
transition_matrix = get_matrix(ENV1 = 0, ENV2 = 0, params, int = 3)
eq = get_eq(pars)

#################################
# MODEL BEHAVIOR
#################################

#data frame
time <- seq(1, 50, by = 0.1)
dm2 <- data.frame(matrix(NA, nrow = length(time), ncol = 4))

#Initial condition
y = c(B = 0.25, T = 0.25, M = 0.25)
dm2[1, ] = c(y, 0.25)
#loop
for(i in 2:length(time)) {
  eq = runsteady(y = y, func = model, parms = pars, times = c(0, i))[[1]]
  eq[4] <- 1 - eq[1] - eq[2] - eq[3]
  dm2[i, ] <- eq
}

#plot
plot(time, dm2[, 1], type = "l", ylim = c(0, max(dm2)), col = 2, ylab = "State proportion")
lines(dm2[, 2], col = 3)
lines(dm2[, 3], col = 4)
lines(dm2[, 4], col = 5)

#################################
# ADD DISTURBANCE
#################################

#data frame with two environmnetal conditions
behavior <- function(envComb1, envComb2) {
	time <- seq(1, 40, by = 1)
	dat <- data.frame(matrix(NA, nrow = length(time), ncol = 6))
	names(dat) <- c("ENV1", "ENV2", "B", "T", "M", "R")
	dat[c(1:20), c(1,2)] <- envComb1 #environmnetal 1
	dat[c(21: dim(dat)[1]), c(1,2)] <- envComb2 #environmnetal 2

	#initial condition
	dat[1, c(3:6)] = c(B = 0.25, T = 0.25, M = 0.25, R = 0.25)

	#behavior
	for(i in 2:length(time)) {
		pars = get_pars(ENV1 = dat[i, 1],ENV2 = dat[i, 2], params, int = 3)
		y = c(B = dat[(i - 1), 3], T = dat[(i - 1), 4], M = dat[(i - 1), 5])
		eq = runsteady(y = y, func = model, parms = pars, times = c(0, i))[[1]]
	  eq[4] <- 1 - eq[1] - eq[2] - eq[3]
	  dat[i, c(3:6)] <- eq
	}

	plot(0, xlim = c(0, length(time)), ylim = c(0, 1), type = 'l', col = 1)
	legend('topleft', c("B", "T", "M", "R"), col = 1:4, lty = 1)

	return(dat)
}

dat <- behavior(envComb1 = c(-.25, -.25), envComb2= c(0.25, 0.25))
lines(dat$B, col = 1)
lines(dat$T, col = 2)
lines(dat$M, col = 3)
lines(dat$R, col = 4)

#################################
# Effect of parameters variation on Eigenvalue
#################################

#TEST 1: original value parameters + 5
#Setting parameters variation (Really ugly)
nn <- 5
By <- 0.5
int <- 2
ab <- seq(from = params[1,], to = params[1,] + nn, By); ab <- 1-(1-logit_reverse(ab))^int
at <- seq(from = params[8,], to = params[8,] + nn, By); at <- 1-(1-logit_reverse(at))^int
bb <- seq(from = params[15,], to = params[15,] + nn, By); bb <- 1-(1-logit_reverse(bb))^int
bt <- seq(from = params[22,], to = params[22,] + nn, By); bt <- 1-(1-logit_reverse(bt))^int
th <- seq(from = params[36,], to = params[36,] + nn, By); th <- 1-(1-logit_reverse(th))^int
tt <- seq(from = params[29,], to = params[29,] + nn, By); tt <- 1-(1-logit_reverse(tt))^int
e <- seq(from = params[43,], to = params[43,] + nn, By); e <- 1-(1-logit_reverse(e))^int
par.var1 <- list(ab, at, bb, bt, th, tt, e)

#running eigenvalue to each parameter
pars = get_pars(ENV1 = 0, ENV2 = 0, params, int = int)
eql <- as.list("NA")
df <- data.frame()
for(k in 1: length(pars)) {
	pars = get_pars(ENV1 = 0, ENV2 = 0, params, int = int)
	for(j in 1: length(ab)) {
		pars[k] = par.var1[[k]][j]
		df[j, 1]	<- par.var1[[k]][j]
		df[j, 2] <- get_eq(pars)$ev
	}
eql[[k]] <- df
}

#plot
par(mfrow = c(3,3), mar = c(2, 2, 1, 1))
for(i in 1:7) {
	plot(eql[[i]], type = "l", xlab = "", ylab = "", ylim = c(-.32,-.02))
	points(eql[[i]])
}

#TEST 2: Fixed value for main parameters from 0 to 1.7 and original value for the other parameters
int <- 2
Fix <- seq(0, 1.7, 0.1)

#running eigenvalue to each parameter
pars = get_pars(ENV1 = 0, ENV2 = 0, params, int = int)
eql <- as.list("NA")
df <- data.frame()
for(k in 1: length(pars)) {
	pars = get_pars(ENV1 = 0, ENV2 = 0, params, int = int)
	for(j in 1: length(Fix)) {
		pars[k] = Fix[j]
		df[j, 1]	<- Fix[j]
		df[j, 2] <- get_eq(pars)$ev
	}
eql[[k]] <- df
}

#plot
par(mfrow = c(3,3), mar = c(2, 2, 1, 1))
for(i in 1:7) {
	plot(eql[[i]], type = "l", xlab = "", ylab = "", ylim = c(-.35,0))
	points(eql[[i]])
}

#TEST 3: Fixed value for all parameters from 0 to 1.7
int <- 2
Fix <- seq(0, 1.7, 0.1)
exPar <- 0.2 #fixed value of all other parameters

#running eigenvalue to each parameter
pars = get_pars(ENV1 = 0, ENV2 = 0, params, int = int)
eql <- as.list("NA")
df <- data.frame()
for(k in 1: length(pars)) {
	pars = get_pars(ENV1 = 0, ENV2 = 0, params, int = int)
	pars[-k] <- exPar
	for(j in 1: length(Fix)) {
		pars[k] = Fix[j]
		df[j, 1]	<- Fix[j]
		df[j, 2] <- get_eq(pars)$ev
	}
eql[[k]] <- df
}

#plot
par(mfrow = c(3,3), mar = c(2, 2, 1, 1))
for(i in 1:7) {
	plot(eql[[i]], type = "l", xlab = "", ylab = "", ylim = c(-.5,0))
	points(eql[[i]])
}
