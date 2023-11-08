
####################################################################################
###                             AUTHOR: DUNG NGOC NGUYEN                         ###
### Exploration of the search space of Gaussian graphical models for paired data ###
####################################################################################


### GREEDY SEARCH ON THE MODEL INCLUSION LATTICE (SECTION 3.4) ###

### A colored graph is composed by components L, E, and E_L.
### In R environment, we encode a colored graph respectively by a list of L.as, E, E.as,
### where L.as is a vector of integers from 1 to p, and E, E_L are two-column matrices with 
### each row corresponding to an edge in the graph, and p is the number of vertices.

### For the model selection procedure on the model inclusion lattice in Section 3.4.1, 
### we have written the function backwardsubmodel() which employs the
### backward elimination stepwise approach for the family of RCON models for paired data 
### with the model inclusion order defined in Section 1.2.4.


### Furthermore, to specify more efficiently the set of neighbors N' described at Step 3.c,
### we wrote the function meet.operation.submodel() which takes the meet operation between 
### the chosen model M* at current step and each (accepted) models in the previous step by. 
### The meet between two models in the model inclusion lattice is identified in Section 1.2.4.


### The procedure can be summarized as follows:

# (1) # From the saturated model M*, we find the set of neighbors of M* by Step 3.a.

# (2) # After testing, we create a vector ii of indices 0 and 1 such that it assigns 0 for 
### rejected models and 1 for accepted models.

# (3) # Update M* by choosing the best model among models having indices 1 in vector ii.

# ... # Going further stages, # (1) # is replaced by specifying the set of neighbors of M*
### through meet.operation.submodel() with inputs M* and models having index 1 in ii.
### Then, we go to # (2) # and # (3) # and repeat the procedure until the stop condition holds.


### backwardsubmodel() returns a list of locally optimal model, number of iterations, 
### number of tested models starting from the iteration.


###========================================================================###
###========================================================================###
###========================================================================###


### LOAD MULTIPLE PACKAGES AT ONCE
load_pkg <- rlang::quos(igraph, gRc, tictoc, foreach)

invisible(lapply(lapply(load_pkg, rlang::quo_name),
                 library,
                 character.only = TRUE
))


### import functions in supplementary_functions.R
source("supplementary_functions.R")


###========================================================================###
###===== MEET OPERATION BETWEEN TWO MODELS USING SUBMODEL RELATION ========###
###===================== IN SECTION 3.4 IN GEHRMANN =======================###
###========================================================================###

### INPUT: 
#         mi: a list of (L.as, E, E.as) for i-th colored graphical model, i = 1, 2
#         p: a (even) numeric (number of vertices)

### OUTPUT: a list of (L.as, E, E.as) of the "meet operation" model, i.e. m1 ^ m2

### [Gehrmann, H. (2011). Lattices of graphical Gaussian models with symmetries]

meet.operation.submodel <- function(m1, m2, p){
  # L.as
  L.as <- sort(intersect(m1$L.as, m2$L.as))
  if(length(L.as) == 0){ L.as <- NULL }
  
  # E
  E1 <- m1$E
  E2 <- m2$E
  E <- intersectMat(E1, E2)$val
  if(!is.null(E)){
    E <- E[order(E[,1],E[,2]),,drop = FALSE]
  }
  
  # E.as
  E.as1 <- m1$E.as
  E.as2 <- m2$E.as
  E.as <- intersectMat(E.as1, E.as2)$val
  if(!is.null(E.as)){
    E.as <- E.as[order(E.as[,1],E.as[,2]),,drop = FALSE]
  }
  
  TE1 <- outEdges(E1, p)$TE
  TE2 <- outEdges(E2, p)$TE
  Es1 <- matdiff(TE1, E.as1)$val
  Es2 <- matdiff(TE2, E.as2)$val
  
  # for m1
  TE1.full <- rbind(TE1, tauMat(TE1, p))
  ET1 <- outEdges(E1, p)$ET
  singleE1 <- matdiff(E1, rbind(ET1, TE1.full))$val
  singleEL1 <- outEdges(singleE1, p)$EL
  singleER1 <- outEdges(singleE1, p)$ER
  ll1 <- intersectMat(singleEL1, Es2)$ind
  rr1 <- intersectMat(tauMat(singleER1, p), Es2)$ind
  E <- matdiff(E, rbind(singleEL1[ll1,, drop = FALSE], singleER1[rr1,, drop = FALSE]))$val
  
  # for m2
  TE2.full <- rbind(TE2, tauMat(TE2, p))
  ET2 <- outEdges(E2, p)$ET
  singleE2 <- matdiff(E2, rbind(ET2, TE2.full))$val
  singleEL2 <- outEdges(singleE2, p)$EL
  singleER2 <- outEdges(singleE2, p)$ER
  ll2 <- intersectMat(singleEL2, Es1)$ind
  rr2 <- intersectMat(tauMat(singleER2, p), Es1)$ind
  E <- matdiff(E, rbind(singleEL2[ll2,, drop = FALSE], singleER2[rr2,, drop = FALSE]))$val
  
  return(list(L.as = L.as, E = E, E.as = E.as))
  
}


###========================================================================###
###============== NEIGHBORHOOD SET OF THE SATURATED MODEL =================###
###========================================================================###

### INPUT: 
#         p: a (even) numeric (number of vertices)
#         m: a list of (L.as, E, E.as) of the saturated model

### OUTPUT: a list of colored graphical models which are the nearest neighbor models
###         of the saturated model.

nei <- function(p, m){
  q <- p/2
  L.as <- m$L.as
  E <- m$E
  E.as <- m$E.as
  
  FT <- fullEdges(p)$FT
  FL <- fullEdges(p)$FL
  
  # vertex
  out <- lapply(1:q, function(i){
    return(list(L.as = setdiff(L.as, i), E = E, E.as = E.as))
  })
  
  # ET
  out <- append(out, lapply(1:q, function(i){
    return(list(L.as = L.as, E = matdiff(E, FT[i,,drop = FALSE])$val, E.as = E.as))
  }))
  
  # FL
  nFL <- nrow(FL)
  temp <- lapply(1:nFL, function(i){
    E.as.f <- rbind(FL[i,], tauMat(FL[i,,drop = FALSE], p))
    mm1 <- list(list(L.as = L.as, E = E, E.as = FL[-i,,drop = FALSE]))
    mm2 <- lapply(1:2, function(j){
      list(L.as = L.as, E = matdiff(E, E.as.f[j,,drop = FALSE])$val,
           E.as = matdiff(E.as, FL[i,,drop = FALSE])$val)
    })
    return(append(mm1, mm2))
  })
  out <- append(out, unlist(temp, recursive = FALSE))
  return(out)
}



###===============================================================================###
###=== THE BACKWARD ELIMINATION STEPWISE PROCEDURE FOR MODEL INCLUSION LATTICE ===###
###===============================================================================###

### INPUT: 
#         df: a data frame
#         itmax: the maximum iterations
#         alpha: the significant level, which is set 0.05 by default
#         parallel: a logic, if it is FALSE (by default) then the binary operator %do%
#         is operated in foreach(), otherwise, %dopar% is applied. This argument might be
#         a numberic, known as positive integer number, which displays the number of cores
#         we would like to run in parallel computation.

### OUTPUT: a list of 
###       model: a list of (L.as, E, E.as) of the optimal model in local search
###       iterations: number of iterations before breaking the procedure
###       no.models: number of tested models which are counted in the iterative manner

backwardsubmodel <- function(df, 
                             itmax,
                             alpha = 0.05,
                             parallel = FALSE){
  # data
  p <- ncol(df)
  q <- p/2
  n <- nrow(df)
  Sig <- cov(df)
  vnames <- c(paste0("L", 1:q, sep=""), paste("R", 1:q, sep=""))
  dimnames(Sig) <- list(vnames, vnames)
  
  FL <- fullEdges(p)$FL
  FV <- fullEdges(p)$FV
  
  # start parallel computations-----------------------------------------------
  if ( !inherits(parallel, "cluster") ) {    # if parallel not already started then
    if ( parallel | is.numeric(parallel) ) {
      parallel <- GA::startParallel(parallel)
      # cluster will be closed here
      on.exit( if (parallel) parallel::stopCluster(attr(parallel, "cluster")) )
    }
  } else {  # if parallel already started
    parallel <- attr(parallel, "cluster")
    # else cluster will be closed in 'emCovGraph' function
  }
  `%DO%` <- if ( is.logical(parallel) ) {
    if ( parallel ) foreach::`%dopar%` else foreach::`%do%`
  } else foreach::`%dopar%`
  #---------------------------------------------------------------------------
  
  # im <- c("intersectMat", "matdiff", "tau", "tauMat", "vLabel", "fullEdges", "outEdges", "rcox.lists")
  # initialization ---------------------------------------------------------
  m0 <- list(L.as = 1:q, E = FV, E.as = FL)
  
  # fitting initial model
  VE.list0 <- rcox.lists(p, m0)
  fit0 <- gRc::rcox(vcc = VE.list0$v.list, ecc = VE.list0$e.list, S = Sig, 
                    n = n, method = "ipms")
  logLik0 <- fit0$fitInfo$logL 
  dof0 <- length(VE.list0$v.list) + length(VE.list0$e.list) 
  # ------------------------------------------------------------------------
  
  ## there are two steps before entering the iterative manner
  
  # step 1  ----------------------------------------------------------------
  S1 <- nei(p, m0)
  pvals.p1 <- foreach::foreach(i = 1:length(S1), .combine=c) %DO% (
    pval.function(p, S1[[i]], logLik0, dof0, Sig, n)
  ) #
  ind.Ap1 <- which(pvals.p1 >= alpha)
  
  ii <- numeric(length(S1))
  ii[ind.Ap1] <- 1
  
  ## choose models to test in step 2
  #
  jj <- ii[(p+1):length(S1)]
  
  #
  s1 <- seq(1, 3*p*(p-2)/4, 3)
  s2 <- seq(3, 3*p*(p-2)/4, 3)
  tt <- foreach::foreach(i = 1:length(s1), .combine=c) %DO% (
    if(all(jj[s1[i]:s2[i]] == 1)){i}
  ) 
  
  # step 2  ----------------------------------------------------------------
  miss.symmE <- FL[tt,, drop = FALSE] 
  if(length(tt) > 0){
    S2 <- foreach::foreach(i = 1:length(tt)) %DO% (
      list(L.as = 1:q, E = matdiff(FV, rbind(miss.symmE[i,], tauMat(miss.symmE[i,, drop = FALSE], p)))$val, 
           E.as = matdiff(FL, miss.symmE[i,,drop = FALSE])$val)
    )
    
    pvals.p2 <- foreach::foreach(i = 1:length(S2), .combine=c) %DO% (
      pval.function(p, S2[[i]], logLik0, dof0, Sig, n)
    ) #
    
    # update S2
    ind.Ap2 <- which(pvals.p2 >= alpha)
    S2 <- S2[ind.Ap2]
    
    # update missing symm EL
    miss.symmE <- miss.symmE[ind.Ap2,, drop = FALSE]
  } else { S2 <- list() }
  # ------------------------------------------------------------------------
  
  
  # choose best model   ----------------------------------------------------
  opt.p <- max(pvals.p1)
  if(opt.p < alpha){
    return(list(model = m0, iterations = 0, no.models = 0))
  } else {
    idxm <- which.max(pvals.p1)
    m <- S1[[idxm]]
  }
  # ------------------------------------------------------------------------
  
  it <- count.models <- 0
  while(it < itmax){
    # from m
    mE <- m$E
    TEm <- outEdges(mE, p)$TE
    Esm <- matdiff(TEm, m$E.as)$val
    ETm <- outEdges(mE, p)$ET
    SE <- matdiff(mE, rbind(rbind(TEm, tauMat(TEm, p)),ETm))$val # single edges
    SL <- outEdges(SE, p)$EL # single edges on the left
    SR <- outEdges(SE, p)$ER # single edges on the right
    ttE1 <- rbind(Esm, rbind(SL, tauMat(SR, p)))
    tt1 <- intersectMat(FL, ttE1)$ind.row
    
    dontTest.ind1 <- foreach::foreach(i = tt1, .combine=c) %DO% ((s1[i]:s2[i]) + p)
    dontTest.ind2 <- c(dontTest.ind1, which(ii == 0))
    
    indToTest.S1 <- setdiff(1:length(S1), unique(c(idxm, dontTest.ind2))) 
    if(length(indToTest.S1) > 0){
      Nei1 <- foreach::foreach(i = indToTest.S1) %DO% (
        meet.operation.submodel(m, S1[[i]], p)
      )
      
      pvals1 <- foreach::foreach(i = 1:length(Nei1), .combine=c) %DO% (
        pval.function(p, Nei1[[i]], logLik0, dof0, Sig, n)
      ) #
      temp.Rp1 <- indToTest.S1[which(pvals1 < alpha)]
      ii[temp.Rp1] <- 0
      #ind.Rp1 <- c(ind.Rp1, temp.Rp1)
      
      # best model's measure
      temp.idxm <- indToTest.S1[which.max(pvals1)]
      opt.p1 <- pvals1[which.max(pvals1)]
    } else {opt.p1 <- 0}
    
    ## search models to test in step 2
    findModel.indS2 <- intersectMat(miss.symmE, ttE1)$ind.row
    if(length(findModel.indS2) > 0){
      Nei2 <- foreach::foreach(i = findModel.indS2) %DO% (
        meet.operation.submodel(m, S2[[i]], p)
      )
      
      pvals2 <- foreach::foreach(i = 1:length(Nei2), .combine=c) %DO% (
        pval.function(p, Nei2[[i]], logLik0, dof0, Sig, n)
      ) #
      
      # best model's pvalue
      opt.p2 <- pvals2[which.max(pvals2)]
      
      # update S2
      ind.Rp2 <- findModel.indS2[which(pvals2 < alpha)]
      if(length(ind.Rp2) > 0){
        # update S2, missing symm EL
        S2 <- S2[-ind.Rp2]
        miss.symmE <- miss.symmE[-ind.Rp2,, drop = FALSE]
      }
      
    } else { opt.p2 <- 0 }
    
    # update missing symm EL from step 1
    jj <- ii[(p+1):length(S1)]
    tt <- foreach::foreach(i = 1:length(s1), .combine=c) %DO% (
      if(all(jj[s1[i]:s2[i]] == 1)){i}
    ) 
    
    #miss.symmE <- intersectMat(miss.symmE, FL[tt,, drop = FALSE])$val 
    diff.S2.ii <- matdiff(miss.symmE, FL[tt,, drop = FALSE])$ind.row
    if(!is.null(diff.S2.ii)) {S2 <- S2[-diff.S2.ii]}
    
    ## count number of models in neighborhood set and iterations
    counting <- length(indToTest.S1) + length(findModel.indS2)
    if(counting > 0){ it <- it + 1 }
    count.models <- count.models + counting
    
    # update m
    if(all(c(opt.p1, opt.p2) < alpha)){
      return(list(model = m, iterations = it, no.models = count.models))
    } else {
      if(opt.p1 > opt.p2){
        m <- Nei1[[which.max(pvals1)]]
        idxm <- c(idxm, temp.idxm)
      } else {
        m <- Nei2[[which.max(pvals2)]]
      }
    }
  }
  return(list(model = m, iterations = it, no.models = count.models))
}

## compile the function
backwardsubmodelc <- compiler::cmpfun(backwardsubmodel)
