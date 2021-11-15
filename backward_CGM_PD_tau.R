#####################################################################################
###                            AUTHOR: DUNG NGOC NGUYEN                           ###
### [Working paper]: Model selection for colored graphical models for paired data ###
#####################################################################################

### A colored graph is composed by components |L in (3.7), the edge set E, |E_L in (3.6)
### in R environment, we translate a colored graph by a corresponding list of L.as, E, E.as


### The model selection procedure in the section 4.1 is described by backwardCGMpd().
### This adapts the idea of backward method in the stepwise approach for the space of RCON models for paired data
### equipped by the partial ordering "tau" defined in the section 3.2. Then, all the operations
### we applied here are induced by the ordering "tau".


### The neighborhood set N' contains models which are submodels of the chosen model 
### in the previous step and they should not be contained by any models in the rejected set R.
### To find this neighborhood set more efficiently, in practise, we do differently
### by applying the meet operation between the chosen model at current step and 
### each (accepted) models in the previous step. We create meet.operation() to do that 
### by making the intersection of sets L.as, E, E.as between two models.

### From the saturated model, we split its neighborhoood set into 2 layers, upper and lower layers.
### layer1() and layer2() are aimed to find the models in upper and lower layers, respectively,
### from the saturated model.

# (1) #
### For the upper layer, we create a vector ii1 to label for each model, after testing, 
### ii1[K] = 0 if the model m_K is rejected, otherwise, ii1[K] = 1.

# (2) #
### For each model m_K in the upper layer satisfying the condition in step 5.a, we will have 
### two corresponding models m_K1 and m_K2 in lower layer. We also create a vector ii2 to label
### for each model in this layer, before testing, we update ii2 by applying the following rule:
### if ii1[K] = 1 then ii2[K1] =  ii2[K2] = 0. Then, we apply the testing for models where (ii2 != 0).
### Update ii2 similarly, ii2[M] = 0 if the model m_M is rejected, otherwise, ii2[M] = 1.

### Going to further stage of the procedure, we might consider models 
### (determined by the fourth point in step 4.a) where # (1) # does not consider them.
### Therefore, to unify the procedure in # (1), (2) #, we apply layer3() to determine models
### which are respectively removed one pair of edges (e, tau(e)) from the full model where
### e != tau(e). For K > p, if ii1[K] = 0 then ii3[K'] = 0 where K'= K-p the corresponding index
### in layer3(). We apply the test for models which have indices in ii3 != 0 and update ii3 
### similarly in # (1), (2) #.

### Starting from the saturated model m*, in the iterative manner, at each stage, 
### the neighboor models in upper layer are the models obtained by taking the meet operation
### between m* and m[ii1 == 1] for m in layer1(). Moreover, if m* satisfies 4th condition in 4.a,
### then m* ^ m[ii3 == 1], for CORRESPONDING m in layer3(). After testing, we update ii1 and ii3,
### and choose the temporary best model, called temp.m1.

### Next, the neighboor models in lower layer of m* are the models obtained by taking
### the meet operation between m* and m[ii2 == 1] for m in layer2() then update ii2 and
### update m* by comparing the best accepted model in layer2() and temp.m1.

### Repeat the procedure until all indices in ii1, ii2, ii3 are 0; or, iterative manner
### exceeds the given maximum number of iterations.


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
###==================== MEET OPERATION IN THEOREM 3.1 =====================###
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
###================== NEIGHBORHOOD SET IN UPPER LAYER =====================###
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
###================== NEIGHBORHOOD SET IN LOWER LAYER =====================###
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
###========================= NEIGHBORHOOD SET 3 ===========================###
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
###= PERFORMANCE OF THE BACKWARD PROCEDURE USING THE PARTIAL ORDERING TAU =###
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
