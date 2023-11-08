######################################################################################
###                              AUTHOR: DUNG NGOC NGUYEN                          ###
#### Exploration of the search space of Gaussian graphical models for paired data ####
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
load("res14.36.tauc.RData")
load("res15.36.tauc.RData")


###=====================================================================================###
###================================ FIGURES 3.10 AND 3.11 ==============================###
###=====================================================================================###

###=====================================================================================###
###             ASYMMETRIC EDGES (FIGURES IN THE FIRST ROW OF 3.10 AND 3.11)            ###
###=====================================================================================###

### SUBJECT 14
outGraph(p, list(L.as =  res14.36.tauc[[1]]$L.as,
                 E = matdiff(outEdges(res14.36.tauc[[1]]$E, 36)$EL, outEdges(res14.36.tauc[[1]]$E, 36)$TE)$val,
                 E.as = NULL), 16)

outGraph(p, list(L.as =  res14.36.tauc[[1]]$L.as,
                 E = matdiff(outEdges(res14.36.tauc[[1]]$E, 36)$ER, tauMat(outEdges(res14.36.tauc[[1]]$E, 36)$TE, p))$val,
                 E.as = NULL), 16)

### SUBJECT 15
outGraph(p, list(L.as = res15.36.tauc[[1]]$L.as, 
                 E = matdiff(outEdges(res15.36.tauc[[1]]$E, 36)$EL, outEdges(res15.36.tauc[[1]]$E, 36)$TE)$val,
                 E.as = NULL), 16)

outGraph(p, list(L.as = res15.36.tauc[[1]]$L.as, 
                 E = matdiff(outEdges(res15.36.tauc[[1]]$E, 36)$ER, tauMat(outEdges(res15.36.tauc[[1]]$E, 36)$TE, p))$val,
                 E.as = NULL), 16)


###=====================================================================================###
###        SYMMETRIC TWIN EDGES (FIGURES IN THE SECOND ROW OF 3.10 AND 3.11)            ###
###=====================================================================================###

### We analyze on the subject 14, and one can apply similarly for subject 15.
p <- 36
q <- p/2

#---------------- within hemispheres

## extract components of colored graph from the selected model
L.as <- res14.36.tauc[[1]]$L.as
E <- res14.36.tauc[[1]]$E
E.as <- res14.36.tauc[[1]]$E.as

## E intersection {(i,j) in VxV: i,j in L}
t1L <- E[E[,1] <= q,,drop=FALSE]
t2L <- t1L[t1L[,2] <= q,,drop=FALSE]

## E intersection {(i,j) in VxV: i,j in R}
t1R <- E[E[,1] <= p & q < E[,1] ,,drop=FALSE]
t2R <- t1R[t1R[,2] <= p & q < t1R[,2],,drop=FALSE]

## E.as intersection {(i,j) in VxV: i,j in L}
t1E.asL <- E.as[E.as[,1] <= q,,drop=FALSE]
t2E.asL <- t1E.asL[t1E.asL[,2] <= q,,drop=FALSE]

## edges within hemispheres 
E.G1 <- rbind(t2L, t2R)
typeG1 <- list(L.as = L.as, E = E.G1, E.as = t2E.asL)

## twin edges
TE.G1 <- outEdges(E.G1, p)$TE

## symmetric twin edges
Es.G1 <- matdiff(TE.G1, t2E.asL)$val
Esf.G1 <- rbind(Es.G1, tauMat(Es.G1,p)) # full edge set of symmetric twin edges

## plotting
par(mfrow=c(1,1), mar=c(0,0,0,0)) 
G1.sym <- list(L.as = L.as, E = Esf.G1, E.as = NULL)
outGraph(p, G1.sym, 16)
#---------------- within hemispheres





#---------------- between hemispheres

## edges between hemispheres
E.G2 <- matdiff(E, E.G1)$val
E.G2 <- E.G2[-which(E.G2[,1] == tau(E.G2[,2], p)),, drop = FALSE] ## minus edges in EI

E.as.G2 <- matdiff(E.as, t2E.asL)$val
typeG2 <- list(L.as = L.as, E = E.G2, E.as = E.as.G2)

## twin edges
TE.G2 <- outEdges(E.G2, p)$TE

## symmetric twin edges
Es.G2 <- matdiff(TE.G2, E.as.G2)$val
Esf.G2 <- rbind(Es.G2, tauMat(Es.G2,p)) # full edge set of symmetric twin edges

## plotting
G2.sym <- list(L.as = L.as, E = Esf.G2, E.as = NULL)
outGraph(p, G2.sym, 16)
#---------------- between hemispheres

