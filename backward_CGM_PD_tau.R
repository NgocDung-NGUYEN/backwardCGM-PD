#########################################################################################
###                              AUTHOR: DUNG NGOC NGUYEN                             ###
### [Doctoral research]: Model selection for colored graphical models for paired data ###
#########################################################################################

### GREEDY SEARCH ON THE TWIN LATTICE (SECTION 3.3) ###

### A colored graph is composed by components L, E, and E_L.
### In R environment, we encode a colored graph respectively by a list of L.as, E, E.as,
### where L.as is a vector of integers from 1 to p, and E, E_L are two-column matrices with 
### each row corresponding to an edge in the graph, and p is the number of vertices.


### For the model selection procedure for twin lattice in Section 3.3.1, we have written
### the function backwardCGMpd() which employs the backward elimination stepwise approach
### for the family of RCON models for paired data with the twin order "tau" defined in Section 2.4.2.
### Then, all the operations are induced by the order "tau".


### Furthermore, to specify more efficiently the sets of neighbors N1' and N2' for the upper layer
### and the lower layer described at Steps 4.c and 5.c, respectively, we wrote the function
### meet.operation() which takes the meet operation between the chosen model M* at current step
### and each (accepted) models in the previous step by the twin order. The meet between two models 
### in the twin lattice is identified in Theorem 2.6.


### The procedure can be summarized as follows:

# (1) # From the saturated model M*, we find the set of neighbors of M* by Steps 4.a and 5.a 
### and split it into 2 layers, upper and lower, by functions layer1() and layer2(), respectively.

# (2) # For the upper layer:
### after testing, we create a vector ii1 of indices 0 and 1 such that it assigns 
### 0 for rejected models and 1 for accepted models.

# (3) # For the lower layer:
### we create a vector ii2 such that for each model having symmetric twin edges in the upper layer
### that is accepted, the corresponding models in the lower layer are assigned to 0, 
### and 1 for others. Then, we only test models in the lower layer having the index 1 and 
### update the vector ii2 by 0 for rejected models.

# (4) # Update M* by choosing the best model among models having indices 1 in vectors ii1 and ii2.

# ... # Going further stages, # (1) # is replaced by specifying the sets of neighbors of M* 
### of the upper and lower layers through meet.operation() with inputs M* and models having 
### index 1 in ii1 and ii2, respectively. Then, we go to # (2) #, # (3) #, and # (4) # 
### and repeat the procedure until the stop condition holds.


#=============================================================================#

### Here, we note that, since the neighborhood of the sarturated models is determined by
### 3 first points in Step 4.a, then for further steps, we wrote layer3() to indicate models
### in the last point of Step 4.a and apply the same manner in # (2), (3) #. layer3() is applied
### when M* has symmetric twin edges.

### backwardCGMpd() returns a list of locally optimal model, number of iterations, 
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
###==================== MEET OPERATION IN THEOREM 2.6 =====================###
###========================================================================###

### INPUT: 
#         mi: a list of (L.as, E, E.as) for i-th colored graphical model, i = 1, 2

### OUTPUT: a list of (L.as, E, E.as) of the "meet operation" model, i.e. m1 ^ m2


# ORDER
meet.operation <- function(m1, m2){
  
  # L.as
  L.as <- sort(intersect(m1$L.as, m2$L.as))
  if(length(L.as) == 0){
    L.as <- NULL
  }
  
  # E
  E <- intersectMat(m1$E, m2$E)$val
  if(!is.null(E)){
    E <- E[order(E[,1],E[,2]),,drop = FALSE]
  }
  
  # E.as
  E.as <- intersectMat(m1$E.as, m2$E.as)$val
  if(!is.null(E.as)){
    E.as <- E.as[order(E.as[,1],E.as[,2]),,drop = FALSE]
  }
  return(list(L.as = L.as, E = E, E.as = E.as))
}


###========================================================================###
###================== SET OF NEIGHBORS IN UPPER LAYER =====================###
###========================================================================###

### INPUT: 
#         p: a (even) numeric (number of vertices)
#         m: a list of (L.as, E, E.as) of the saturated model

### OUTPUT: a list of colored graphical models which are the nearest neighbor models
###         of the saturated model in the upper layer 

layer1 <- function(p, m){
  L.as <- m$L.as
  E <- m$E
  E.as <- m$E.as
  
  ET <- outEdges(E, p)$ET
  
  q <- p/2
  
  #
  out <- lapply(1:q, function(i){
    return(list(L.as = setdiff(L.as, i), E = E, E.as = E.as))
  })
  
  #
  out <- append(out, lapply(1:nrow(ET), function(i){
    list(L.as = L.as, E = matdiff(E, ET[i,,drop = FALSE])$val, E.as = E.as)
  }))
  
  #
  out <- append(out, lapply(1:nrow(E.as), function(i){
    list(L.as = L.as, E = E, E.as = E.as[-i,,drop = FALSE])
  }))
  
  return(out)
}


###========================================================================###
###================== SET OF NEIGHBORS IN LOWER LAYER =====================###
###========================================================================###

### INPUT: 
#         p: a (even) numeric (number of vertices)
#         m: a list of (L.as, E, E.as) of the saturated model
#         ind.R: a vector of indices indicating rejected models in the upper layer
#                which have indices greater than p. ind.R contains extracted indices
#                after minus them by p.

### OUTPUT: a list of colored graphical models which are the nearest neighbor models
###         of the saturated model in the lower layer 

layer2 <- function(p, m, ind.R){
  L.as <- m$L.as
  E <- m$E
  E.as <- m$E.as
  
  out <- list()
  idx <- ind.R[ind.R > p] - p
  
  
  if(length(idx) != 0){
    temp <- lapply(idx, function(i){
      E.as.f <- rbind(E.as[i,], tauMat(E.as[i,,drop = FALSE], p))
      lapply(1:2, function(j){
        list(L.as = L.as, E = matdiff(E, E.as.f[j,,drop = FALSE])$val,
             E.as = E.as[-i,, drop = FALSE])
      })
    })
    out <- append(out, unlist(temp, recursive = FALSE))
  }
  
  return(out)
}


###========================================================================###
###============ SET OF MODELS IN THE LAST POINT OF STEP 4.a ===============###
###========================================================================###

### INPUT: 
#         p: a (even) numeric (number of vertices)
#         m: a list of (L.as, E, E.as) of the saturated model
#         ind.A: a vector of indices indicating accepted models in the upper layer
#                which have indices greater than p. ind.A contains extracted indices
#                after minus them by p.

### OUTPUT: a list of colored graphical models which are removed respectively each pair
###         of edges (e, tau(e)) from the saturated model.

layer3 <- function(p, m, ind.A){
  L.as <- m$L.as
  E <- m$E
  E.as <- m$E.as
  
  out <- list()
  idx <- ind.A[ind.A > p] - p
  
  if(length(idx) != 0){
    out <- append(out, lapply(idx, function(i){
      Es.fm <- rbind(E.as[i,], tauMat(E.as[i,, drop = FALSE], p))
      return(list(L.as = L.as, E = matdiff(E, Es.fm)$val, E.as = E.as[-i,, drop = FALSE]))
    }))
  }
  
  return(out)
}


###========================================================================###
###== THE BACKWARD ELIMINATION STEPWISE PROCEDURE FOR THE TWIN LATTICE  ===###
###========================================================================###

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
 

backwardCGMpd <- function(df,
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
  
  # test models in layer 1  ------------------------------------------------
  L1 <- layer1(p, m0)
  pvals.p1 <- foreach::foreach(i = 1:length(L1), .combine=c) %DO% (
    pval.function(p, L1[[i]], logLik0, dof0, Sig, n)
  ) #
  ind.Ap1 <- which(pvals.p1 >= alpha)
  ind.Rp1 <- which(pvals.p1 < alpha)
  
  ii1 <- numeric(length(L1))
  ii1[ind.Ap1] <- 1
  
  # 
  idx.TE.from.A1 <- ind.Ap1[ind.Ap1 > p] - p
  TE.from.A1 <- FL[idx.TE.from.A1,, drop = FALSE]
  # ------------------------------------------------------------------------
  
  # test models in layer 2  ------------------------------------------------
  L2 <- layer2(p, m0, ind.Rp1)
  if(length(L2) > 0){
    pvals.p2 <- foreach::foreach(i = 1:length(L2), .combine=c) %DO% (
      pval.function(p, L2[[i]], logLik0, dof0, Sig, n)
    ) #
    ind.Ap2 <- which(pvals.p2 >= alpha)
    
    ii2 <- numeric(length(L2))
    ii2[ind.Ap2] <- 1 
  } else{ii2 <-  pvals.p2 <- 0 }
  # ------------------------------------------------------------------------
  
  # test models in layer 3  ------------------------------------------------
  L3 <- layer3(p, m0, ind.Ap1)
  if(length(L3) > 0){
    pvals.p3 <- foreach::foreach(i = 1:length(L3), .combine=c) %DO% (
      pval.function(p, L3[[i]], logLik0, dof0, Sig, n)
    ) #
    ind.Ap3 <- which(pvals.p3 >= alpha)
    ii3 <- numeric(length(L3))
    ii3[ind.Ap3] <- 1
  } else { ii3 <- 0}
  # ------------------------------------------------------------------------
  
  
  # choose best model from layer 1 and 2 -----------------------------------
  opt.p1 <- max(pvals.p1)
  opt.p2 <- max(pvals.p2)
  
  if(max(opt.p1, opt.p2) < alpha){
    return(list(model = m0, iterations = 0, no.models = 0))
  } else {
    if(opt.p1 >= opt.p2){
      idxm.A1 <- which.max(pvals.p1)
      next.idxm.A2 <- idxm.A2 <- idxm.A3 <- c()
      m <- L1[[idxm.A1]]
    } else {
      idxm.A2 <- which.max(pvals.p2)
      next.idxm.A2 <- ifelse(idxm.A2 %% 2 == 0, idxm.A2 - 1, idxm.A2 + 1)
      idxm.A1 <- idxm.A3 <- c()
      m <- L2[[idxm.A2]]
    }
  }
  # ------------------------------------------------------------------------
  
  it <- count.models <- 0
  while(it < itmax){
    
    ### layer 1
    #
    temp.A1 <- setdiff(1:length(L1), c(idxm.A1, which(ii1 == 0)))
    if(length(temp.A1) > 0){
      Nei1.model <- foreach::foreach(i = temp.A1) %DO% ( meet.operation(m, L1[[i]])) #
      
      pvals1 <- foreach::foreach(i = 1:length(Nei1.model), .combine=c) %DO% (
        pval.function(p, Nei1.model[[i]], logLik0, dof0, Sig, n)
      ) #
      
      ii1[temp.A1[which(pvals1 < alpha)]] <- 0
      opt.p11 <- max(pvals1)
    } else { opt.p11 <- 0 } # ok
    
    #
    TE.m <- outEdges(m$E, p)$TE
    Es.m <- matdiff(TE.m, m$E.as)$val
    idx.intersect <- intersectMat(TE.from.A1, Es.m)$ind.row
    temp.A3 <- setdiff(idx.intersect, c(idxm.A3, which(ii3 == 0)))
    if(length(temp.A3) > 0){
      Nei3.model <- foreach::foreach(i = temp.A3) %DO% (meet.operation(m, L3[[i]])) #
      
      pvals3 <- foreach::foreach(i = 1:length(Nei3.model), .combine=c) %DO% (
        pval.function(p, Nei3.model[[i]], logLik0, dof0, Sig, n)
      ) #
      
      ii3[temp.A3[which(pvals3 < alpha)]] <- 0
      opt.p13 <- max(pvals3)
    } else { opt.p13 <- 0 }
    
    ### choose the best model in layer 1
    if(max(opt.p11, opt.p13) == 0){
      opt.p1 <- 0
    } else {
      if(opt.p11 > opt.p13){
        opt.i <- which.max(pvals1)
        tt1 <- temp.A1[opt.i]
        opt.p1 <- opt.p11
        temp.m1 <- Nei1.model[[opt.i]]
        temp.idxm <- 1
      } else {
        opt.i <- which.max(pvals3)
        tt1 <- temp.A3[opt.i]
        opt.p1 <- opt.p13
        temp.m1 <- Nei3.model[[opt.i]]
        temp.idxm <- 3
      }
    }
    
    ### layer 2
    dontTest.ind <- c(idxm.A2, next.idxm.A2)
    temp.A2 <- setdiff(1:length(L2), unique(c(dontTest.ind, which(ii2 == 0))))
    if(length(temp.A2) > 0 & all(temp.A2 == 0) == FALSE){
      Nei2.model <- foreach::foreach(i = temp.A2) %DO% (meet.operation(m, L2[[i]])) #
      
      pvals2 <- foreach::foreach(i = 1:length(Nei2.model), .combine=c) %DO% (
        pval.function(p, Nei2.model[[i]], logLik0, dof0, Sig, n)
      ) #
      ii2[temp.A2[which(pvals2 < alpha)]] <- 0
      opt.p2 <- max(pvals2)
    } else { opt.p2 <- 0 }
    
    ## count number of models in neighborhood set and iterations
    counting <- length(temp.A1) + length(temp.A3) + length(temp.A2)
    if(counting > 0){ it <- it + 1 }
    count.models <- count.models + counting
    print(it)
    
    ## update m (compare two best models in layer 1 and layer 2)
    if(max(opt.p1, opt.p2) < alpha){
      return(list(model = m, iterations = it, no.models = count.models))
    } else {
      if(opt.p1 >= opt.p2){
        m <- temp.m1
        if(temp.idxm == 1){idxm.A1 <- c(idxm.A1, tt1) } else {idxm.A3 <- c(idxm.A3, tt1)}
      } else {
        opt.i <- which.max(pvals2)
        tt2 <- temp.A2[opt.i]
        m <- Nei2.model[[opt.i]]
        idxm.A2 <- c(idxm.A2, tt2)
        next.idxm.A2 <- c(next.idxm.A2, ifelse(tt2 %% 2 == 0, tt2 - 1, tt2 + 1))
      }
    }
  }
  return(list(model = m, iterations = it, no.models = count.models))
}

## compile the function
backwardCGMpdc <- compiler::cmpfun(backwardCGMpd)
