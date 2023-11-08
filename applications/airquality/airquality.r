#####################################################################################
###                            AUTHOR: DUNG NGOC NGUYEN                           ###
#### Exploration of the search space of Gaussian graphical models for paired data ###
#####################################################################################

### LOAD MULTIPLE PACKAGES AT ONCE
load_pkg <- rlang::quos(igraph, gRc, tictoc, foreach, pdglasso)

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
load("airdataAR1.RData")
n <- nrow(data.res)
p <- ncol(data.res)

#------------------------------------------------------------------------#
### sample variance 
S <- var(data.res)
plot(diag(S), pch=as.character(c(1:6, 1:6)))

### fit model using pd-glasso
mod.S <- pdRCON.fit(S, n, gamma.eBIC = 0)
GS <- pdColG.get(mod.S$model)


#------------------------------------------------------------------------#
# sample correlation
P <- cov2cor(S)
plot(diag(P), pch=as.character(c(1:6, 1:6)))

### fit model using pd-glasso
mod.P <- pdRCON.fit(P, n, gamma.eBIC = 0)
GP <- pdColG.get(mod.P$model)


#------------------------------------------------------------------------#
### compute p-values of models based on the log-likelihood ratio test
p <- ncol(data.res)
q <- p/2
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

### pvalues of GS and GP models

#
GS1 <- GS$pdColG
colnames(GS1) <- rownames(GS1) <- NULL
GS.mod <- formation.matrix(p, NULL, GS1)
pval.function(p, GS.mod, logLik0, dof0, Sig, n, Kstart = NULL) 

# 
GP1 <- GP$pdColG
colnames(GP1) <- rownames(GP1) <- NULL
GP.mod <- formation.matrix(p, NULL, GP1)
pval.function(p, GP.mod, logLik0, dof0, Sig, n, Kstart = NULL) 


#------------------------------------------------------------------------#
### airquality in backward procedure
air.backward <- backwardCGMpd(data.res, itmax = 1000, alpha = 0.05, parallel = 3)

### pvalue
pval.function(p, air.backward[[1]] , logLik0, dof0, Sig, n, Kstart = NULL)


