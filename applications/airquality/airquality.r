#####################################################################################
###                            AUTHOR: DUNG NGOC NGUYEN                           ###
#### Exploration of the search space of Gaussian graphical models for paired data ###
#####################################################################################

### LOAD MULTIPLE PACKAGES AT ONCE
load_pkg <- rlang::quos(igraph, gRc, tictoc, foreach)

invisible(lapply(lapply(load_pkg, rlang::quo_name),
                 library,
                 character.only = TRUE
))

### import sources
# "../" points to the parent folder
### import functions
source("../supplementary_functions.R")
source("../backward_CGM_PD_tau.R")
source("../airsupp.R")


### load dataset
load("airdata.RData")
nv <- dim(data)[2] # number of variables


### check the autocorrelation of the univariate time series 
for(i in 1:nv){
  acf(data[, i], main=names(data[i]))
  readline(prompt = "Press [enter] for the next plot: ")
}

### check the autocorrelation of the residuals obtained by 
### applying an AR(1) model to the univariate series
for(i in 1:nv){
  acf(arima(data[, i], order=c(1,0,0))$residuals, main=names(data[i]))
  readline(prompt = "Press [enter] for the next plot: ")
}

### construct a dataframe with residual of AR(1) model
data.res <- matrix(NA, nrow=373, ncol=nv)
for(i in 1:nv){
  data.res[, i] <- arima(data[, i], order=c(1,0,0))$residuals
}

data.res  <- as.data.frame(data.res) 
names(data.res) <- names(data)

### otherwise, we can load residual of AR(1) model of airdata directly
# load("airdataAR1.RData")

#------------------------------------------------------------------------#
library(pdglasso)

#------------------------------------------------------------------------#
### sample variance 
S <- var(data.res)
n <- dim(data.res)[1]
plot(diag(S), pch=as.character(c(1:6, 1:6)))

### fit model using pd-glasso
mod.S <- pdRCON.fit(S, n, gamma.eBIC = 0)
pdRCON.check(mod.S$model)

#
GS <- pdColG.get(mod.S$model)


#------------------------------------------------------------------------#
# sample correlation
P <- cov2cor(S)
n <- dim(data.res)[1]
plot(diag(P), pch=as.character(c(1:6, 1:6)))

### fit model using pd-glasso
mod.P <- pdRCON.fit(P, n, gamma.eBIC = 0)
pdRCON.check(mod.P$model)

#
GP <- pdColG.get(mod.P$model)


#------------------------------------------------------------------------#
### compute p-values of models based on the log-likelihood ratio test
p <- ncol(data.res)
q <- p/2
n <- nrow(data.res)
Sig <- var(data.res)
vnames <- c(paste0("L", 1:q, sep=""), paste("R", 1:q, sep=""))
dimnames(Sig) <- list(vnames, vnames)

FL <- fullEdges(p)$FL
FV <- fullEdges(p)$FV

### saturated model
m0 <- list(L.as = 1:q, E = FV, E.as = FL)

### fit the saturated model
VE.list0 <- rcox.lists(p, m0)
fit0 <- gRc::rcox(vcc = VE.list0$v.list, ecc = VE.list0$e.list, S = Sig, 
                  n = n, method = "ipms")

### log-likelihood value
logLik0 <- fit0$fitInfo$logL 

### number of degrees
dof0 <- length(VE.list0$v.list) + length(VE.list0$e.list) 
Kstart <- NULL

### pvalues of GS and GP models

#
GS1 <- GS$pdColG
colnames(GS1) <- rownames(GS1) <- NULL
GS.mod <- formation.matrix(p, NULL, GS1)
pval.function(p, GS.mod, logLik0, dof0, Sig, n, Kstart = Kstart) 
outGraph(p, GS.mod, 13)

# 
GP1 <- GP$pdColG
colnames(GP1) <- rownames(GP1) <- NULL
GP.mod <- formation.matrix(p, NULL, GP1)
pval.function(p, GP.mod, logLik0, dof0, Sig, n, Kstart = Kstart) 
outGraph(p, GP.mod, 13)


#------------------------------------------------------------------------#
### airquality in backward procedure
air.backward <- backwardCGMpd(data.res, itmax = 1000, alpha = 0.05, parallel = 3)


### pvalue
pval.function(p, air.backward[[1]] , logLik0, dof0, Sig, n, Kstart = Kstart)
outGraph(p,  air.backward[[1]], 13)

#------------------------------------------------------------------------#
### Matrix representations
par(mfrow = c(1, 1))

# GS
plotG(GS$pdColG, p, NULL)

# GP
plotG(GP$pdColG, p, NULL) 

# air.backward
mat <- formation.matrix(p, air.backward[[1]], NULL)
plotG(mat, p, NULL) 
