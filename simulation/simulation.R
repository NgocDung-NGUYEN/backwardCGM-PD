
######################################################################################
###                              AUTHOR: DUNG NGOC NGUYEN                          ###
### [Ph.D. research]: Model selection for colored graphical models for paired data ###
######################################################################################

### LOAD MULTIPLE PACKAGES AT ONCE
load_pkg <- rlang::quos(igraph, gRc, tictoc, foreach)

invisible(lapply(lapply(load_pkg, rlang::quo_name),
                 library,
                 character.only = TRUE
))

### import sources
# "../" points to the parent folder
source("../supplementary_functions.R")
source("../backward_CGM_PD_tau.R")
source("../backward_CGM_PD_submod.R")


########################################################################################
############################## SIMULATION Section 3.5.1 ################################
########################################################################################

################################################################################
### p |||        |||             A              |||            B             ###
###==========================================================================###
### 8 |||  edges  ||| 1 sym, 1 non-sym (pairs)  ||| 3 sym, 1 non-sym (pairs) ###
###   |||vertices ||| 1 sym, 3 non-sym (pairs)  ||| 3 sym, 1 non-sym (pairs) ###
###--------------------------------------------------------------------------###
### 12 ||| edges  ||| 1 sym, 4 non-sym (pairs)  ||| 6 sym, 3 non-sym (pairs) ###
###    |||vertices||| 2 sym, 4 non-sym (pairs)  ||| 4 sym, 2 non-sym (pairs) ###
###--------------------------------------------------------------------------###
### 16 ||| edges  ||| 2 sym, 7 non-sym (pairs)  ||| 10 sym, 5 non-sym (pairs)###
###    |||vertices||| 2 sym, 6 non-sym (pairs)  ||| 6 sym, 2 non-sym (pairs) ###
###--------------------------------------------------------------------------###
### 20 ||| edges  ||| 3 sym, 11 non-sym (pairs) ||| 16 sym, 8 non-sym (pairs)###
###    |||vertices||| 2 sym, 8 non-sym (pairs)  ||| 8 sym, 2 non-sym (pairs) ###
################################################################################

### Note that, var.11.8 refers to a variable in scenario A with the graph of 8 vertices,
### and var.22.8 refers to a variable in scenario B with the graph of 8 vertices.

########################################################################################

## number of observations
n <- 100

## rho in covariance matrix
rho <- 0.5


##################
### SCENARIO A ###
##################


#############
### p = 8 ###
#############
p11.8 <- 8
q11.8 <- p11.8/2

FV <- fullEdges(p11.8)$FV
FL <- fullEdges(p11.8)$FL
FI <- fullEdges(p11.8)$FI


## equicovariance matrix S
eq.S11.8 <- matrix(rho, ncol = p11.8, nrow = p11.8) 
diag(eq.S11.8) <- 1
vlabs11.8 <- c(paste0("L", 1:q11.8, sep=""), paste("R", 1:q11.8, sep=""))
dimnames(eq.S11.8) <- list(vlabs11.8, vlabs11.8)

## E.as
set.seed(9)
pairsE11.8 <- FL[sample(1:nrow(FL), 2),] # 2 pairs
pairsE11.8.full <- rbind(pairsE11.8, tauMat(pairsE11.8, p11.8))
E.as11.8 <- pairsE11.8[1,, drop = FALSE] # 1 pairs of non-sym

## E
restFV11.8 <- matdiff(FV, rbind(pairsE11.8.full, FI))$val
set.seed(1)
restE11.8 <- restFV11.8[sample(1:nrow(restFV11.8), 1),, drop = FALSE]

E11.8 <- rbind(pairsE11.8.full, restE11.8)

## L
set.seed(1)
v11.8 <- sort(sample(1:q11.8, 3))

## colored graph G
g11.8 <- list(L.as = v11.8, E = E11.8, E.as = E.as11.8)
outGraph(p11.8, g11.8, 8)

## fit the concentration matrix from eq. matrix and G
cc11.8 <- rcox.lists(p11.8, g11.8)
gfit11.8 <- gRc::rcox(vcc = cc11.8$v.list, ecc = cc11.8$e.list, 
                      S = eq.S11.8, n = n, method = "ipms")
Theta11.8 <- gfit11.8$fitInfo$K


## check positive definite
Sigma11.8 <- solve(Theta11.8)
all(eigen(Sigma11.8)$value > 0) 

## generate data
simdf11.8 <- lapply(1:20, function(i) MASS::mvrnorm(n = n, mu = rep(0, p11.8), Sigma = Sigma11.8))


##############
### p = 12 ###
##############
p11.12 <- 12
q11.12 <- p11.12/2

FV <- fullEdges(p11.12)$FV
FL <- fullEdges(p11.12)$FL
FI <- fullEdges(p11.12)$FI


## equicovariance matrix S
eq.S11.12 <- matrix(rho, ncol = p11.12, nrow = p11.12) 
diag(eq.S11.12) <- 1
vlabs11.12 <- c(paste0("L", 1:q11.12, sep=""), paste("R", 1:q11.12, sep=""))
dimnames(eq.S11.12) <- list(vlabs11.12, vlabs11.12)

## EI
set.seed(1)
EI11.12 <- FI[sample(1:nrow(FI), 1),] #1 

## E.as
set.seed(2)
pairsE11.12 <- FL[sample(1:nrow(FL), 5),] # 5 pairs
pairsE11.12.full <- rbind(pairsE11.12, tauMat(pairsE11.12, p11.12))
E.as11.12 <- pairsE11.12[sample(1:nrow(pairsE11.12), 4),] # 4 pairs of non-sym

## E
restFV11.12 <- matdiff(FV, rbind(pairsE11.12.full, FI))$val
set.seed(1)
restE11.12 <- restFV11.12[sample(1:nrow(restFV11.12), 2),]

E11.12 <- rbind(pairsE11.12.full, restE11.12)

## L
set.seed(1)
v11.12 <- sort(sample(1:q11.12, 4))

## colored graph G
g11.12 <- list(L.as = v11.12, E = E11.12, E.as = E.as11.12)
outGraph(p11.12, g11.12, 23)

## fit the concentration matrix from eq. matrix and G
cc11.12 <- rcox.lists(p11.12, g11.12)
gfit11.12 <- gRc::rcox(vcc = cc11.12$v.list, ecc = cc11.12$e.list, 
                       S = eq.S11.12, n = n, method = "ipms")
Theta11.12 <- gfit11.12$fitInfo$K


## check positive definite
Sigma11.12 <- solve(Theta11.12)
all(eigen(Sigma11.12)$value > 0) 

## generate data
simdf11.12 <- lapply(1:20, function(i) MASS::mvrnorm(n = n, 
                                                     mu = rep(0, p11.12), 
                                                     Sigma = Sigma11.12))

###=====================================================================================###
###=====================================================================================###
###=====================================================================================###


##############
### p = 16 ###
##############
p11.16 <- 16
q11.16 <- p11.16/2

FV <- fullEdges(p11.16)$FV
FL <- fullEdges(p11.16)$FL
FI <- fullEdges(p11.16)$FI


## equicovariance matrix S
eq.S11.16 <- matrix(rho, ncol = p11.16, nrow = p11.16) 
diag(eq.S11.16) <- 1
vlabs11.16 <- c(paste0("L", 1:q11.16, sep=""), paste("R", 1:q11.16, sep=""))
dimnames(eq.S11.16) <- list(vlabs11.16, vlabs11.16)

## EI
set.seed(1)
EI11.16 <- FI[sample(1:nrow(FI), 1),, drop = FALSE] #2 EI

## E.as
set.seed(2)
pairsE11.16 <- FL[sample(1:nrow(FL), 9),] # 9 pairs
pairsE11.16.full <- rbind(pairsE11.16, tauMat(pairsE11.16, p11.16))
set.seed(2)
E.as11.16 <- pairsE11.16[sample(1:nrow(pairsE11.16), 7),] # 7 pairs of non-sym

## E
restFV11.16 <- matdiff(FV, rbind(pairsE11.16.full, FI))$val
set.seed(1)
restE11.16 <- restFV11.16[sample(1:nrow(restFV11.16), 3),]

E11.16 <- rbind(rbind(EI11.16, pairsE11.16.full), restE11.16)

## L
set.seed(1)
v11.16 <- sort(sample(1:q11.16, 6))

## colored graph G
g11.16 <- list(L.as = v11.16, E = E11.16, E.as = E.as11.16)
outGraph(p11.16, g11.16, 5)

## fit the concentration matrix from eq. matrix and G
cc11.16 <- rcox.lists(p11.16, g11.16)
gfit11.16 <- gRc::rcox(vcc = cc11.16$v.list, ecc = cc11.16$e.list, 
                       S = eq.S11.16, n = n, method = "ipms")
Theta11.16 <- gfit11.16$fitInfo$K


## check positive definite
Sigma11.16 <- solve(Theta11.16)
all(eigen(Sigma11.16)$value > 0) 

## generate data
simdf11.16 <- lapply(1:20, function(i) MASS::mvrnorm(n = n, 
                                                     mu = rep(0, p11.16), 
                                                     Sigma = Sigma11.16))

###=====================================================================================###
###=====================================================================================###
###=====================================================================================###


##############
### p = 20 ###
##############
p11.20 <- 20
q11.20 <- p11.20/2

FV <- fullEdges(p11.20)$FV
FL <- fullEdges(p11.20)$FL
FI <- fullEdges(p11.20)$FI


## equicovariance matrix S
eq.S11.20 <- matrix(rho, ncol = p11.20, nrow = p11.20) 
diag(eq.S11.20) <- 1
vlabs11.20 <- c(paste0("L", 1:q11.20, sep=""), paste("R", 1:q11.20, sep=""))
dimnames(eq.S11.20) <- list(vlabs11.20, vlabs11.20)

## EI
set.seed(1)
EI11.20 <- FI[sample(1:nrow(FI), 2),] #2 EI

## E.as
set.seed(2)
pairsE11.20 <- FL[sample(1:nrow(FL), 14),] # 14 pairs
pairsE11.20.full <- rbind(pairsE11.20, tauMat(pairsE11.20, p11.20))

set.seed(2)
E.as11.20 <- pairsE11.20[sample(1:nrow(pairsE11.20), 11),] # 11 pairs of non-sym

## E
restFV11.20 <- matdiff(FV, rbind(pairsE11.20.full, FI))$val
set.seed(4)
restE11.20 <- restFV11.20[sample(1:nrow(restFV11.20), 4),]
outEdges(restE11.20, 20)$TE

E11.20 <- rbind(rbind(EI11.20, pairsE11.20.full), restE11.20)

## L
set.seed(1)
v11.20 <- sort(sample(1:q11.20, 8))

## colored graph G
g11.20 <- list(L.as = v11.20, E = E11.20, E.as = E.as11.20)
outGraph(p11.20, g11.20, 13)

## fit the concentration matrix from eq. matrix and G
cc11.20 <- rcox.lists(p11.20, g11.20)
gfit11.20 <- gRc::rcox(vcc = cc11.20$v.list, ecc = cc11.20$e.list, 
                       S = eq.S11.20, n = n, method = "ipms")
Theta11.20 <- gfit11.20$fitInfo$K


## check positive definite
Sigma11.20 <- solve(Theta11.20)
all(eigen(Sigma11.20)$value > 0) 

## generate data
simdf11.20 <- lapply(1:20, function(i) MASS::mvrnorm(n = n, 
                                                     mu = rep(0, p11.20), 
                                                     Sigma = Sigma11.20))


###=====================================================================================###
###========================== RUN SIMMULATION FOR SCENARIO A ===========================###
###=====================================================================================###


#############
### p = 8 ###
#############

## get the results
out.11.8.tau <- out.11.8.submod <- vector(mode = "list", length = 20)
n.steps.11.8.tau <- n.steps.11.8.submod <- c()
n.models.11.8.tau <- n.models.11.8.submod <- c()
runtime.11.8.tau <- runtime.11.8.submod <- c()

tic() # 
for (i in 1:20) {
  print(i)
  
  #
  ptm <- proc.time()
  f1 <- backwardCGMpd1(simdf11.8[[i]], itmax = 500, alpha = 0.05, parallel = 3)
  t1 <- proc.time() - ptm
  print(t1)
  runtime.11.8.tau <- c(runtime.11.8.tau, as.numeric(t1)[3])
  out.11.8.tau[[i]] <- f1[[1]]
  n.steps.11.8.tau <- c(n.steps.11.8.tau, f1[[2]])
  n.models.11.8.tau <- c(n.models.11.8.tau, f1[[3]])
  
  #
  ptm <- proc.time()
  f2 <- backwardsubmodel(simdf11.8[[i]], itmax = 500, alpha = 0.05, parallel = 3)
  t2 <- proc.time() - ptm
  print(t2)
  runtime.11.8.submod <- c(runtime.11.8.submod, as.numeric(t2)[3])
  out.11.8.submod[[i]] <- f2[[1]]
  n.steps.11.8.submod <- c(n.steps.11.8.submod, f2[[2]])
  n.models.11.8.submod <- c(n.models.11.8.submod, f2[[3]])
  
}
toc()

###=====================================================================================###

##############
### p = 12 ###
##############

## get the results
out.11.12.tau <- out.11.12.submod <- vector(mode = "list", length = 20)
n.steps.11.12.tau <- n.steps.11.12.submod <- c()
n.models.11.12.tau <- n.models.11.12.submod <- c()
runtime.11.12.tau <- runtime.11.12.submod <- c()

tic() # 2548.896 sec elapsed
for (i in 1:20) {
  print(i)
  
  #
  ptm <- proc.time()
  f1 <- backwardCGMpd1(simdf11.12[[i]], itmax = 500, alpha = 0.05, parallel = 3)
  t1 <- proc.time() - ptm
  print(t1)
  runtime.11.12.tau <- c(runtime.11.12.tau, as.numeric(t1)[3])
  out.11.12.tau[[i]] <- f1[[1]]
  n.steps.11.12.tau <- c(n.steps.11.12.tau, f1[[2]])
  n.models.11.12.tau <- c(n.models.11.12.tau, f1[[3]])
  
  #
  ptm <- proc.time()
  f2 <- backwardsubmodel(simdf11.12[[i]], itmax = 500, alpha = 0.05, parallel = 3)
  t2 <- proc.time() - ptm
  print(t2)
  runtime.11.12.submod <- c(runtime.11.12.submod, as.numeric(t2)[3])
  out.11.12.submod[[i]] <- f2[[1]]
  n.steps.11.12.submod <- c(n.steps.11.12.submod, f2[[2]])
  n.models.11.12.submod <- c(n.models.11.12.submod, f2[[3]])
  
}
toc()

###=====================================================================================###

##############
### p = 16 ###
##############

## get the results
out.11.16.tau <- out.11.16.submod <- vector(mode = "list", length = 20)
n.steps.11.16.tau <- n.steps.11.16.submod <- c()
n.models.11.16.tau <- n.models.11.16.submod <- c()
runtime.11.16.tau <- runtime.11.16.submod <- c()

tic() # 12423.132 sec elapsed
for (i in 1:20) {
  print(i)
  
  #
  ptm <- proc.time()
  f1 <- backwardCGMpd1(simdf11.16[[i]], itmax = 500, alpha = 0.05, parallel = 3)
  t1 <- proc.time() - ptm
  print(t1)
  runtime.11.16.tau <- c(runtime.11.16.tau, as.numeric(t1)[3])
  out.11.16.tau[[i]] <- f1[[1]]
  n.steps.11.16.tau <- c(n.steps.11.16.tau, f1[[2]])
  n.models.11.16.tau <- c(n.models.11.16.tau, f1[[3]])
  
  #
  ptm <- proc.time()
  f2 <- backwardsubmodel(simdf11.16[[i]], itmax = 500, alpha = 0.05, parallel = 3)
  t2 <- proc.time() - ptm
  print(t2)
  runtime.11.16.submod <- c(runtime.11.16.submod, as.numeric(t2)[3])
  out.11.16.submod[[i]] <- f2[[1]]
  n.steps.11.16.submod <- c(n.steps.11.16.submod, f2[[2]])
  n.models.11.16.submod <- c(n.models.11.16.submod, f2[[3]])
  
}
toc()

###=====================================================================================###

##############
### p = 20 ###
##############

## get the results
out.11.20.tau <- out.11.20.submod <- vector(mode = "list", length = 20)
n.steps.11.20.tau <- n.steps.11.20.submod <- c()
n.models.11.20.tau <- n.models.11.20.submod <- c()
runtime.11.20.tau <- runtime.11.20.submod <- c()

tic() # 2:24 (22006.478 sec elapsed)
for (i in 19:20) {
  print(i)
  
  #
  ptm <- proc.time()
  f1 <- backwardCGMpd1(simdf11.20[[i]], itmax = 500, alpha = 0.05, parallel = 3)
  t1 <- proc.time() - ptm
  print(t1)
  runtime.11.20.tau <- c(runtime.11.20.tau, as.numeric(t1)[3])
  out.11.20.tau[[i]] <- f1[[1]]
  n.steps.11.20.tau <- c(n.steps.11.20.tau, f1[[2]])
  n.models.11.20.tau <- c(n.models.11.20.tau, f1[[3]])
  
  #
  ptm <- proc.time()
  f2 <- backwardsubmodel(simdf11.20[[i]], itmax = 500, alpha = 0.05, parallel = 3)
  t2 <- proc.time() - ptm
  print(t2)
  runtime.11.20.submod <- c(runtime.11.20.submod, as.numeric(t2)[3])
  out.11.20.submod[[i]] <- f2[[1]]
  n.steps.11.20.submod <- c(n.steps.11.20.submod, f2[[2]])
  n.models.11.20.submod <- c(n.models.11.20.submod, f2[[3]])
  
}
toc()


###=====================================================================================###
###=====================================================================================###
###=====================================================================================###

### PLOT by ggplot() for averaged elapsed time

runtime.df.tau <- data.frame(c(p11.8, p11.12, p11.16, p11.20), 
                             c(mean(runtime.11.8.tau), mean(runtime.11.12.tau),
                               mean(runtime.11.16.tau),mean(runtime.11.20.tau)))
runtime.df.submod <- data.frame(c(p11.8, p11.12, p11.16, p11.20), 
                                c(mean(runtime.11.8.submod), mean(runtime.11.12.submod),
                                  mean(runtime.11.16.submod),mean(runtime.11.20.submod)))
colnames(runtime.df.tau) <- colnames(runtime.df.submod) <- c("p", "seconds")

par(mfrow=c(1,1), mar=c(0,0,0,0))
ggplot() + 
  geom_point(data=runtime.df.tau,aes(x=p, y=seconds, color="tau"), size=2.5, shape=19)  + 
  geom_line(data=runtime.df.tau,aes(x=p, y=seconds), color="red", size=1)  + 
  geom_point(data=runtime.df.submod, aes(x=p, y=seconds, color="submod"), size=2.5, shape=19) +
  geom_line(data=runtime.df.submod,aes(x=p, y=seconds), color="blue", size=1)  + 
  scale_color_manual("", breaks=c("tau", "submod"), values = c("red","blue"))+
  scale_x_continuous(breaks = c(8, 12, 16, 20),
                     labels = c("8", "12", "16", "20"))+
  ylab("average seconds") + xlab("number of vertices p")+
  ggtitle("Running time for scenario A")

###=====================================================================================###

### PLOT by ggplot() for averaged numbers of fitted models in iterative manner
nmodels.df.tau <- data.frame(c(p11.8, p11.12, p11.16, p11.20), 
                             c(mean(n.models.11.8.tau), mean(n.models.11.12.tau),
                               mean(n.models.11.16.tau),mean(n.models.11.20.tau)))
nmodels.df.submod <- data.frame(c(p11.8, p11.12, p11.16, p11.20), 
                                c(mean(n.models.11.8.submod), mean(n.models.11.12.submod),
                                  mean(n.models.11.16.submod),mean(n.models.11.20.submod)))
colnames(nmodels.df.tau) <- colnames(nmodels.df.submod) <- c("p", "models")

par(mfrow=c(1,1), mar=c(0,0,0,0))
ggplot() + 
  geom_point(data=nmodels.df.tau,aes(x=p, y=models, color="tau"), size=2.5, shape=19)  + 
  geom_line(data=nmodels.df.tau,aes(x=p, y=models), color="red", size=1)  + 
  geom_point(data=nmodels.df.submod, aes(x=p, y=models, color="submod"), size=2.5, shape=19) +
  geom_line(data=nmodels.df.submod,aes(x=p, y=models), color="blue", size=1)  + 
  scale_color_manual("", breaks=c("tau", "submod"), values = c("red","blue"))+
  scale_x_continuous(breaks = c(8, 12, 16, 20),
                     labels = c("8", "12", "16", "20"))+
  ylab("number of fitted models") + xlab("number of vertices p")+
  ggtitle("Scenario A")








###=====================================================================================###
###=====================================================================================###
###=====================================================================================###



##################
### SCENARIO B ###
##################


#############
### p = 8 ###
#############
p22.8 <- 8
q22.8 <- p22.8/2

FV <- fullEdges(p22.8)$FV
FL <- fullEdges(p22.8)$FL
FI <- fullEdges(p22.8)$FI


## equicovariance matrix S
eq.S22.8 <- matrix(rho, ncol = p22.8, nrow = p22.8) 
diag(eq.S22.8) <- 1
vlabs22.8 <- c(paste0("L", 1:q22.8, sep=""), paste("R", 1:q22.8, sep=""))
dimnames(eq.S22.8) <- list(vlabs22.8, vlabs22.8)

## E.as
set.seed(4)
pairsE22.8 <- FL[sample(1:nrow(FL), 4),] # 4 pairs
pairsE22.8.full <- rbind(pairsE22.8, tauMat(pairsE22.8, p22.8))
set.seed(1)
E.as22.8 <- pairsE22.8[sample(1:nrow(pairsE22.8), 1),, drop = FALSE] # 1 pairs of non-sym

## E
restFV22.8 <- matdiff(FV, rbind(pairsE22.8.full, FI))$val
set.seed(2)
restE22.8 <- restFV22.8[sample(1:nrow(restFV22.8), 2),]

E22.8 <- rbind(pairsE22.8.full, restE22.8)

## L
set.seed(1)
v22.8 <- sort(sample(1:q22.8, 1))

## colored graph G
g22.8 <- list(L.as = v22.8, E = E22.8, E.as = E.as22.8)
outGraph(p22.8, g22.8, 8)

### fit the covariance matrix from eq. matrix S and G
cc22.8 <- rcox.lists(p22.8, g22.8)
gfit22.8 <- gRc::rcox(vcc = cc22.8$v.list, ecc = cc22.8$e.list, 
                      S = eq.S22.8, n = n, method = "ipms")
Theta22.8 <- gfit22.8$fitInfo$K


## check positive defnite
Sigma22.8 <- solve(Theta22.8)
all(eigen(Sigma22.8)$value > 0) 

## generate data
simdf22.8 <- lapply(1:20, function(i) MASS::mvrnorm(n = n, mu = rep(0, p22.8), Sigma = Sigma22.8))


###=====================================================================================###
###=====================================================================================###
###=====================================================================================###

##############
### p = 12 ###
##############
p22.12 <- 12
q22.12 <- p22.12/2

FV <- fullEdges(p22.12)$FV
FL <- fullEdges(p22.12)$FL
FI <- fullEdges(p22.12)$FI


## equicovariance matrix S
eq.S22.12 <- matrix(rho, ncol = p22.12, nrow = p22.12) 
diag(eq.S22.12) <- 1
vlabs22.12 <- c(paste0("L", 1:q22.12, sep=""), paste("R", 1:q22.12, sep=""))
dimnames(eq.S22.12) <- list(vlabs22.12, vlabs22.12)

## EI
set.seed(1)
EI22.12 <- FI[sample(1:nrow(FI), 1),, drop = FALSE] #1 EI

## E.as
set.seed(2)
pairsE22.12 <- FL[sample(1:nrow(FL), 9),] # 9 pairs
pairsE22.12.full <- rbind(pairsE22.12, tauMat(pairsE22.12, p22.12))
E.as22.12 <- pairsE22.12[sample(1:nrow(pairsE22.12), 3),] # 3 pairs of non-sym


## E
restFV22.12 <- matdiff(FV, rbind(pairsE22.12.full, FI))$val
set.seed(3)
restE22.12 <- restFV22.12[sample(1:nrow(restFV22.12), 4),, drop = FALSE]

E22.12 <- rbind(rbind(EI22.12, pairsE22.12.full), restE22.12)

## L
set.seed(1)
v22.12 <- sort(sample(1:q22.12, 2))

## colored graph G
g22.12 <- list(L.as = v22.12, E = E22.12, E.as = E.as22.12)
outGraph(p22.12, g22.12, 23)

## fit the covariance matrix from eq. matrix S and G
cc22.12 <- rcox.lists(p22.12, g22.12)
gfit22.12 <- gRc::rcox(vcc = cc22.12$v.list, ecc = cc22.12$e.list, 
                       S = eq.S22.12, n = n, method = "ipms")
Theta22.12 <- gfit22.12$fitInfo$K


## check positive definite
Sigma22.12 <- solve(Theta22.12)
all(eigen(Sigma22.12)$value > 0) 

## generate data
simdf22.12 <- lapply(1:20, function(i) MASS::mvrnorm(n = n, mu = rep(0, p22.12), Sigma = Sigma22.12))

###=====================================================================================###
###=====================================================================================###
###=====================================================================================###


##############
### p = 16 ###
##############
p22.16 <- 16
q22.16 <- p22.16/2

FV <- fullEdges(p22.16)$FV
FL <- fullEdges(p22.16)$FL
FI <- fullEdges(p22.16)$FI


## equicovariance matrix S
eq.S22.16 <- matrix(rho, ncol = p22.16, nrow = p22.16) 
diag(eq.S22.16) <- 1
vlabs22.16 <- c(paste0("L", 1:q22.16, sep=""), paste("R", 1:q22.16, sep=""))
dimnames(eq.S22.16) <- list(vlabs22.16, vlabs22.16)

## EI
set.seed(1)
EI22.16 <- FI[sample(1:nrow(FI), 3),] #3 EI

## E.as
set.seed(2)
pairsE22.16 <- FL[sample(1:nrow(FL), 15),] # 15 pairs
pairsE22.16.full <- rbind(pairsE22.16, tauMat(pairsE22.16, p22.16))
set.seed(2)
E.as22.16 <- pairsE22.16[sample(1:nrow(pairsE22.16), 5),] # 5 pairs of non-sym

## E
restFV22.16 <- matdiff(FV, rbind(pairsE22.16.full, FI))$val
set.seed(1)
restE22.16 <- restFV22.16[sample(1:nrow(restFV22.16), 9),]

E22.16 <- rbind(rbind(EI22.16, pairsE22.16.full), restE22.16)

## L
set.seed(1)
v22.16 <- sort(sample(1:q22.16, 2))

## colored graph G
g22.16 <- list(L.as = v22.16, E = E22.16, E.as = E.as22.16)
outGraph(p22.16, g22.16, 4)

## fit the covariance matrix from eq. matrix S and G
cc22.16 <- rcox.lists(p22.16, g22.16)
gfit22.16 <- gRc::rcox(vcc = cc22.16$v.list, ecc = cc22.16$e.list, 
                       S = eq.S22.16, n = n, method = "ipms")
Theta22.16 <- gfit22.16$fitInfo$K


## check positive definite
Sigma22.16 <- solve(Theta22.16)
all(eigen(Sigma22.16)$value > 0) 

## generate data
simdf22.16 <- lapply(1:20, function(i) MASS::mvrnorm(n = n, mu = rep(0, p22.16), Sigma = Sigma22.16))


###=====================================================================================###
###=====================================================================================###
###=====================================================================================###


##############
### p = 20 ###
##############
p22.20 <- 20
q22.20 <- p22.20/2

FV <- fullEdges(p22.20)$FV
FL <- fullEdges(p22.20)$FL
FI <- fullEdges(p22.20)$FI


## equicovariance matrix S
eq.S22.20 <- matrix(rho, ncol = p22.20, nrow = p22.20) 
diag(eq.S22.20) <- 1
vlabs22.20 <- c(paste0("L", 1:q22.20, sep=""), paste("R", 1:q22.20, sep=""))
dimnames(eq.S22.20) <- list(vlabs22.20, vlabs22.20)

## EI
set.seed(1)
EI22.20 <- FI[sample(1:nrow(FI), 6),] #6 EI

## E.as
set.seed(2)
pairsE22.20 <- FL[sample(1:nrow(FL), 24),] # 24 pairs
pairsE22.20.full <- rbind(pairsE22.20, tauMat(pairsE22.20, p22.20))

set.seed(2)
E.as22.20 <- pairsE22.20[sample(1:nrow(pairsE22.20), 8),] # 8 pairs of non-sym

## E
restFV22.20 <- matdiff(FV, rbind(pairsE22.20.full, FI))$val
set.seed(1)
restE22.20 <- restFV22.20[sample(1:nrow(restFV22.20), 12),]
outEdges(restE22.20, 20)$TE

E22.20 <- rbind(rbind(EI22.20, pairsE22.20.full), restE22.20)

## L
set.seed(1)
v22.20 <- sort(sample(1:q22.20, 2))

## colored graph G
g22.20 <- list(L.as = v22.20, E = E22.20, E.as = E.as22.20)
outGraph(p22.20, g22.20, 21)

## fit the covariance matrix from eq. matrix S and G
cc22.20 <- rcox.lists(p22.20, g22.20)
gfit22.20 <- gRc::rcox(vcc = cc22.20$v.list, ecc = cc22.20$e.list, 
                       S = eq.S22.20, n = n, method = "ipms")
Theta22.20 <- gfit22.20$fitInfo$K


## check positive defininte
Sigma22.20 <- solve(Theta22.20)
all(eigen(Sigma22.20)$value > 0) 

## generate data
simdf22.20 <- lapply(1:20, function(i) MASS::mvrnorm(n = n, mu = rep(0, p22.20), Sigma = Sigma22.20))



###=====================================================================================###
###=========================== RUN SIMMULATION FOR SCENARIO B ==========================###
###=====================================================================================###
#############
### p = 8 ###
#############

## get the results
out.22.8.tau <- out.22.8.submod <- vector(mode = "list", length = 20)
n.steps.22.8.tau <- n.steps.22.8.submod <- c()
n.models.22.8.tau <- n.models.22.8.submod <- c()
runtime.22.8.tau <- runtime.22.8.submod <- c()

tic() # 352.578 sec elapsed
for (i in 1:20) {
  print(i)
  
  #
  ptm <- proc.time()
  f1 <- backwardCGMpd1(simdf22.8[[i]], itmax = 500, alpha = 0.05, parallel = 3)
  t1 <- proc.time() - ptm
  print(t1)
  runtime.22.8.tau <- c(runtime.22.8.tau, as.numeric(t1)[3])
  out.22.8.tau[[i]] <- f1[[1]]
  n.steps.22.8.tau <- c(n.steps.22.8.tau, f1[[2]])
  n.models.22.8.tau <- c(n.models.22.8.tau, f1[[3]])
  
  #
  ptm <- proc.time()
  f2 <- backwardsubmodel(simdf22.8[[i]], itmax = 500, alpha = 0.05, parallel = 3)
  t2 <- proc.time() - ptm
  print(t2)
  runtime.22.8.submod <- c(runtime.22.8.submod, as.numeric(t2)[3])
  out.22.8.submod[[i]] <- f2[[1]]
  n.steps.22.8.submod <- c(n.steps.22.8.submod, f2[[2]])
  n.models.22.8.submod <- c(n.models.22.8.submod, f2[[3]])
  
}
toc()

###=====================================================================================###

##############
### p = 12 ###
##############

## get the results
out.22.12.tau <- out.22.12.submod <- vector(mode = "list", length = 20)
n.steps.22.12.tau <- n.steps.22.12.submod <- c()
n.models.22.12.tau <- n.models.22.12.submod <- c()
runtime.22.12.tau <- runtime.22.12.submod <- c()

tic() # 2428.333 sec elapsed
for (i in 1:20) {
  print(i)
  
  #
  ptm <- proc.time()
  f1 <- backwardCGMpd1(simdf22.12[[i]], itmax = 500, alpha = 0.05, parallel = 3)
  t1 <- proc.time() - ptm
  print(t1)
  runtime.22.12.tau <- c(runtime.22.12.tau, as.numeric(t1)[3])
  out.22.12.tau[[i]] <- f1[[1]]
  n.steps.22.12.tau <- c(n.steps.22.12.tau, f1[[2]])
  n.models.22.12.tau <- c(n.models.22.12.tau, f1[[3]])
  
  #
  ptm <- proc.time()
  f2 <- backwardsubmodel(simdf22.12[[i]], itmax = 500, alpha = 0.05, parallel = 3)
  t2 <- proc.time() - ptm
  print(t2)
  runtime.22.12.submod <- c(runtime.22.12.submod, as.numeric(t2)[3])
  out.22.12.submod[[i]] <- f2[[1]]
  n.steps.22.12.submod <- c(n.steps.22.12.submod, f2[[2]])
  n.models.22.12.submod <- c(n.models.22.12.submod, f2[[3]])
  
}
toc()

###=====================================================================================###

##############
### p = 16 ###
##############

## get the results
out.22.16.tau <- out.22.16.submod <- vector(mode = "list", length = 20)
n.steps.22.16.tau <- n.steps.22.16.submod <- c()
n.models.22.16.tau <- n.models.22.16.submod <- c()
runtime.22.16.tau <- runtime.22.16.submod <- c()

tic() # 12380.526 sec elapsed
for (i in 1:20) {
  print(i)
  
  #
  ptm <- proc.time()
  f1 <- backwardCGMpd1(simdf22.16[[i]], itmax = 500, alpha = 0.05, parallel = 3)
  t1 <- proc.time() - ptm
  print(t1)
  runtime.22.16.tau <- c(runtime.22.16.tau, as.numeric(t1)[3])
  out.22.16.tau[[i]] <- f1[[1]]
  n.steps.22.16.tau <- c(n.steps.22.16.tau, f1[[2]])
  n.models.22.16.tau <- c(n.models.22.16.tau, f1[[3]])
  
  #
  ptm <- proc.time()
  f2 <- backwardsubmodel(simdf22.16[[i]], itmax = 500, alpha = 0.05, parallel = 3)
  t2 <- proc.time() - ptm
  print(t2)
  runtime.22.16.submod <- c(runtime.22.16.submod, as.numeric(t2)[3])
  out.22.16.submod[[i]] <- f2[[1]]
  n.steps.22.16.submod <- c(n.steps.22.16.submod, f2[[2]])
  n.models.22.16.submod <- c(n.models.22.16.submod, f2[[3]])
  
}
toc()

###=====================================================================================###

##############
### p = 20 ###
##############

## get the results
out.22.20.tau <- out.22.20.submod <- vector(mode = "list", length = 20)
n.steps.22.20.tau <- n.steps.22.20.submod <- c()
n.models.22.20.tau <- n.models.22.20.submod <- c()
runtime.22.20.tau <- runtime.22.20.submod <- c()

tic() # 00:12 (53561.885 sec elapsed)
for (i in 1:20) {
  print(i)
  
  #
  ptm <- proc.time()
  f1 <- backwardCGMpd1(simdf22.20[[i]], itmax = 500, alpha = 0.05, parallel = 3)
  t1 <- proc.time() - ptm
  print(t1)
  runtime.22.20.tau <- c(runtime.22.20.tau, as.numeric(t1)[3])
  out.22.20.tau[[i]] <- f1[[1]]
  n.steps.22.20.tau <- c(n.steps.22.20.tau, f1[[2]])
  n.models.22.20.tau <- c(n.models.22.20.tau, f1[[3]])
  
  #
  ptm <- proc.time()
  f2 <- backwardsubmodel(simdf22.20[[i]], itmax = 500, alpha = 0.05, parallel = 3)
  t2 <- proc.time() - ptm
  print(t2)
  runtime.22.20.submod <- c(runtime.22.20.submod, as.numeric(t2)[3])
  out.22.20.submod[[i]] <- f2[[1]]
  n.steps.22.20.submod <- c(n.steps.22.20.submod, f2[[2]])
  n.models.22.20.submod <- c(n.models.22.20.submod, f2[[3]])
  
}
toc()

###=====================================================================================###
###=====================================================================================###
###=====================================================================================###


## PLOT by ggplot for averaged elapsed time

runtime.df.tau22 <- data.frame(c(p22.8, p22.12, p22.16, p22.20), 
                               c(mean(runtime.22.8.tau), mean(runtime.22.12.tau),
                                 mean(runtime.22.16.tau),mean(runtime.22.20.tau)))
runtime.df.submod22 <- data.frame(c(p22.8, p22.12, p22.16, p22.20), 
                                  c(mean(runtime.22.8.submod), mean(runtime.22.12.submod),
                                    mean(runtime.22.16.submod),mean(runtime.22.20.submod)))
colnames(runtime.df.tau22) <- colnames(runtime.df.submod22) <- c("p", "seconds")

ggplot() + 
  geom_point(data=runtime.df.tau22,aes(x=p, y=seconds, color="tau"), size=2.5, shape=19)  + 
  geom_line(data=runtime.df.tau22,aes(x=p, y=seconds), color="red", size=1)  + 
  geom_point(data=runtime.df.submod22, aes(x=p, y=seconds, color="submod"), size=2.5, shape=19) +
  geom_line(data=runtime.df.submod22,aes(x=p, y=seconds), color="blue", size=1)  + 
  scale_color_manual("", breaks=c("tau", "submod"), values = c("red","blue"))+
  scale_x_continuous(breaks = c(8, 12, 16, 20),
                     labels = c("8", "12", "16", "20"))+
  ylab("average seconds") + xlab("number of vertices p")+
  ggtitle("Running time for scenario B")


###=====================================================================================###

## PLOT by ggplot for averaged number of fitted models

nmodels.df.tau22 <- data.frame(c(p22.8, p22.12, p22.16, p22.20), 
                               c(mean(n.models.22.8.tau), mean(n.models.22.12.tau),
                                 mean(n.models.22.16.tau),mean(n.models.22.20.tau)))
nmodels.df.submod22 <- data.frame(c(p22.8, p22.12, p22.16, p22.20), 
                                  c(mean(n.models.22.8.submod), mean(n.models.22.12.submod),
                                    mean(n.models.22.16.submod),mean(n.models.22.20.submod)))
colnames(nmodels.df.tau22) <- colnames(nmodels.df.submod22) <- c("p", "models")

par(mfrow=c(1,1), mar=c(0,0,0,0))
ggplot() + 
  geom_point(data=nmodels.df.tau22,aes(x=p, y=models, color="tau"), size=2.5, shape=19)  + 
  geom_line(data=nmodels.df.tau22,aes(x=p, y=models), color="red", size=1)  + 
  geom_point(data=nmodels.df.submod22, aes(x=p, y=models, color="submod"), size=2.5, shape=19) +
  geom_line(data=nmodels.df.submod22,aes(x=p, y=models), color="blue", size=1)  + 
  scale_color_manual("", breaks=c("tau", "submod"), values = c("red","blue"))+
  scale_x_continuous(breaks = c(8, 12, 16, 20),
                     labels = c("8", "12", "16", "20"))+
  ylab("number of fitted models") + xlab("number of vertices p")+
  ggtitle("Scenario B")


###=====================================================================================###
###================================ RESULTS FOR TABLE 3.2 ==============================###
###=====================================================================================###

## Here we take an example for Scenario B with p = 8

## Users can use the saved version in .RData format of the output if one want to 
## avoid running the code again. This version can be found in simulation-results file.
load("out.22.8.tau.RData")

temp <- out.22.8.tau # here we can change for another case
eA <- E22.8

##
eTP <- sapply(sapply(temp,`[`,2), function(A) {nrow(intersectMat(A, eA)$val)})
nedges <- sapply(sapply(temp,`[`,2), function(A) nrow(A))

##
ePPV <- mean(eTP/nedges) 
ePPV 

##
no.edges <- nrow(eA)
eTPR <- mean(eTP/no.edges) 
eTPR 

## 
FV <- fullEdges(p22.8)$FV

##
miss.E <- lapply(sapply(temp ,`[`,2), function(A) matdiff(FV, A)$val )
miss.E.true <- matdiff(FV, eA)$val
eTN <- sapply(miss.E, function(A) {nrow(intersectMat(A, miss.E.true)$val)})

##
tot.edges <- nrow(FV)
eTNR <- mean(eTN/(tot.edges - no.edges))
eTNR 

##
p <- p22.8
e.asA <- E.as22.8
TE <- outEdges(eA, p)$TE
Es.true <- matdiff(TE, e.asA)$val
TE.list <- lapply(sapply(temp,`[`,2), function(A) outEdges(A, p)$TE)
Es.list <- mapply( function(A,B) matdiff(A,B)$val, TE.list, sapply(temp,`[`,3))
sTP <- sapply(Es.list, function(A) {nrow(intersectMat(A, Es.true)$val)})
nsymm <- sapply(Es.list, function(A) nrow(A))
sP <- nrow(Es.true)

## if sP[i] == NULL then sP[i] = 0
for(i in 1:20){
  if(is.null(sTP[[i]])) {sTP[[i]] <- 0}
}

## if nsymm[i] == NULL then nsymm = 0
for(i in 1:20){
  if(is.null(nsymm[[i]])) {nsymm[[i]] <- 0}
}

## compute sPPV when sP has some elements that are NULL
sPPV <- mean(unlist(sTP)/nsymm)
sPPV #

## compute sTPR when sP has some elements that are NULL
sTPR <- mean(unlist(sTP)/sP)
sTPR # 0.876

## compute sPPV when there is no element that are NULL in sP
sPPV <- mean(sTP/nsymm)
sPPV 

## compute sTPR when there is no element that are NULL in sP
sTPR <- mean(sTP/sP)
sTPR 

##
miss.TE.true <- outEdges(miss.E.true, p)$TE
miss.TE <- lapply(miss.E, function(A) outEdges(A, p)$TE)

##
sTN <- sapply(miss.TE, function(A) nrow(intersectMat(A, miss.TE.true)$val))
sN <- dim(miss.TE.true)[1]

##
sTNR <- mean(sTN/sN)
sTNR







