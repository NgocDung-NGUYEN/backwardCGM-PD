######################################################################################
###                              AUTHOR: DUNG NGOC NGUYEN                          ###
#### Exploration of the search space of Gaussian graphical models for paired data ####
######################################################################################


### A colored graph is composed by components |L in (3.7), the edge set E, |E_L in (3.6)
### in R environment, we translate a colored graph by a corresponding list of L.as, E, E.as


### Let E(G) be a edge set of the graph G, E(G) is performed by a two-column matrix.
### Then, intersectMat() displays the intersection between E(G1) and E(G2);
###       matdiff() displays E(G1)\E(G2).


### tau() shows the twin correspondence in 3.1 and 
### tauMat() is the twin correspondence for edge set.


### In paired data problem, the graph G has an even number of vertices p and 
### we set q = p/2. Then, V = {1,..., p} is partitioned into two disjoint subsets L and R.
### Given a set of vertices V, and p, 
### vLabel() returns a vector of labels of vertices ("L1",..., "Lq") union ("R1",..., "Rq")
### where ("L1",..., "Lq"), ("R1",..., "Rq") are sets of vertices in L, R, respectively.


### fullEdges() returns sets FV, FL, FR, FT, given p, shown in (3.2)-(3.4).


### outEdges() returns sets EL, ER, ET, TE (pairs of twin edges), 
### given edge set E and p, shown in (3.5).


### rcox.lists() returns list of v.list (list of labels of vertex coloring) 
### and e.list (list of labels of edge coloring).


### pval.function() returns the pvalue from the chi-square of the likelihood-ratio testing 
### of a model relative to the full model m0


### FOR VISULIZATION OF A MODEL BY COLORED GRAPH:
### layoutSym() puts coordinates of vertices such that they are positioned symmetrically
### outGraph() display a colored graph from input of elements (L.as, E, E.as)
### where L.as = |L in (3.7), E = edge set, E.as = |E_L in (3.6)


###========================================================================###
###========================================================================###
###========================================================================###


### LOAD MULTIPLE PACKAGES AT ONCE
load_pkg <- rlang::quos(igraph, gRc, tictoc, foreach)

invisible(lapply(lapply(load_pkg, rlang::quo_name),
                 library,
                 character.only = TRUE
))


###========================================================================###
###=================== INTERSECTION OF TWO MATRICES =======================###
###========================================================================###


### INPUT: A, B are two two-column matrices of positive integer numbers 

### OUTPUT: 
#         ind.row: a vector of row indices of A, which contain row entries appearing A and B
#         val: a (sub)matrix of A and B where entries are extracted from ind.row

intersectMat <- function(A,B){
  if(is.null(A) | is.null(B)){
    return(list(ind.row = NULL, val = NULL))
  } else {
    tA <- A[order(A[,1],A[,2]),,drop = FALSE] # temp A
    tB <- B[order(B[,1],B[,2]),,drop = FALSE] # temp B
    
    if(identical(tA,tB)) {
      return(list(ind.row = 1:nrow(A), val = A))
    } else {
      dfA <- data.frame(A)
      dfB <- data.frame(B)
      colsToUse <- colnames(dfA)
      ind <- na.omit(match(do.call("paste", dfB[,colsToUse]), do.call("paste", dfA[,colsToUse])))
      
      if(length(ind) == 0){
        ind <- val  <- NULL
      } else {
        ind <- sort(ind)
        val <- A[ind,, drop = FALSE]
      }
      return(list(ind.row = ind, val = val))
    }
  }
}


###========================================================================###
###==================== DIFFERENCE BETWEEN TWO MATRICES ===================###
###========================================================================###


### INPUT: A, B are two two-column matrices of positive integer numbers 

### OUTPUT:
#         ind.row: a vector of row indices of A, which contain row entries in A and not in B
#         val: a (sub)matrix of A where entries are extracted from ind.row

matdiff <- function(A,B){
  if(is.null(A)) {
    return(list(ind.row = NULL, val = NULL))
  } else if (is.null(B)) {
    return(list(ind.row = 1:nrow(A), val = A))
  } else if(!is.null(A) & !is.null(B)){
    
    temp <- intersectMat(A,B)$ind.row
    ind <- setdiff(1:nrow(A), temp)
    val <- A[ind,, drop = FALSE] 
    if(is.matrix(val) & nrow(val) == 0){
      ind <- val <- NULL
    }
    return(list(ind.row = ind, val = val))
  }
}



###========================================================================###
###=========================== TAU FUNCTION ===============================###
###========================================================================###

### INPUT: 
#         v: a vector or a (natural) numeric (set of vertices)
#         p: a (even) numeric (number of vertices)

### OUTPUT: a vector of tau(v)


tau <- function(v, p){
  if(is.null(v)){
    return(NULL)
  } else {
    q <- p/2
    out <- sapply(v, function(x) {
      if(x %in% 1:q){x + q} else if(x %in% (q+1):p){x - q}
    } )
    return(out)
  }
}


###========================================================================###
###========================== TAU OF A MATRIX =============================###
###========================================================================###

### INPUT: 
#         A: a two-column matrix of positive integer numbers (set of edges)
#         p: a (even) numeric (number of vertices)

### OUTPUT: a two-column matrix of tau(A)

tauMat <- function(A, p){
  if(is.null(A)){
    return(NULL)
  } else {
    if(nrow(A) == 0){ return(NULL)}
    out <- c()
    for (k in 1:nrow(A)) {
      if(tau(A[k,1], p) < tau(A[k,2], p)){
        out <- rbind(out, c(tau(A[k,1], p), tau(A[k,2], p)) )
      } else {
        out <- rbind(out, c(tau(A[k,2], p), tau(A[k,1], p)) )
      }
    }
    return(out)
  }
}


###========================================================================###
###=============== LABELS OF TWO SETS OF HOMOLOGOUS VERICES ===============###
###========================================================================###

### INPUT: 
#         v: a vector or a numeric, interpreted as a positive integer, (set of vertices)
#         p: a (even) numeric (number of vertices)

### OUTPUT: a vector of characters ("L1",..., "Lq") and ("R1",..., "Rq") where q = p/2

vLabel <- function(v, p){
  if(is.null(v)){ return(NULL) } 
  else {
    q <- p/2
    v.lab <- sapply(v, function(x) {
      if(x %in% 1:q){ paste("L", x, sep = "") } 
      else if(x %in% (q+1):p){ paste("R", x-q, sep = "")}
    } ) 
    return(v.lab)
  }
}


###========================================================================###
###=========================== FV, FL, FR, FT =============================###
###========================================================================###

### INPUT: 
#         p: a (even) numeric (number of vertices)

### OUTPUT: a list of
#         FV: a two-column matrix of all possible edges
#         FL: a two-column matrix of all possible edges on the left
#         FR: a two-column matrix of all possible edges on the right
#         FT: a two-column matrix of all possible edges linking (i, tau(i)), i = 1,..,q


fullEdges <- function(p){
  if(p > 0 & (p %% 2) == 0){
    
    # FV
    FV <- expand.grid(x = 1:p, y = 1:p) #Cartesian product (1:p)x(1:p)
    FV <- subset(FV, x < y) #all of possible edges
    
    # FT
    FT <- subset(FV, y == tau(x,p))
    FT <- data.matrix(FT)
    
    # FL
    FL <- subset(FV, (x < tau(1,p) & y < tau(1,p)) | x < tau(y,p))
    FL <- data.matrix(FL)
    
    # FR
    FR <- subset(FV, (x >= tau(1,p) & y >= tau(1,p)) | (x > tau(y,p) & tau(y,p) > 0))
    FR <- data.matrix(FR)
    
    # return all of outputs into matrices
    FV <- data.matrix(FV)
    dimnames(FV) <- dimnames(FT) <- dimnames(FL) <- dimnames(FR) <- NULL
    
    # return all of entries as numeric 
    # (if we don't, it will return integer and make error when using identical())
    FV <- apply(FV, 2, as.numeric)
    FL <- apply(FL, 2, as.numeric)
    FR <- apply(FR, 2, as.numeric)
    FT <- apply(FT, 2, as.numeric)
    
    # order all entries by the first and second column
    FV <- FV[order(FV[,1],FV[,2]),,drop = FALSE]
    FL <- FL[order(FL[,1],FL[,2]),,drop = FALSE]
    FR <- FR[order(FR[,1],FR[,2]),,drop = FALSE]
    FT <- FT[order(FT[,1],FT[,2]),,drop = FALSE]
    
  } else {
    FV <- FT <- FL <- FR <- NULL
  }
  return(list(FV = FV, FT = FT, FL = FL, FR = FR))
}


###========================================================================###
###================= EL, ER, ET, TE (pairs of twin edges) =================###
###========================================================================###

### INPUT: 
#         E: a two-column matrix of edges
#         p: a (even) numeric (number of vertices)

### OUTPUT: a list of
#         EL: a two-column matrix of edges on the left
#         ER: a two-column matrix of edges on the right
#         ET: a two-column matrix of edges linking (i, tau(i)), for some i = 1,..,q
#         TE: a two-column matrix of pairs of edges ( in EL intersection tau(ER) )


# ORDER (from function fullEdges())
outEdges <- function(E, p){
  if(!is.null(E)){
    all.edges <- fullEdges(p)
    EL <- intersectMat(all.edges$FL, E)$val
    ER <- intersectMat(all.edges$FR, E)$val
    ET <- intersectMat(all.edges$FT, E)$val
    twin.E <- intersectMat(EL, tauMat(ER, p))$val 
  } else {
    EL <- ER <- ET <- twin.E <- NULL
  }
  return(list(ET = ET, EL = EL, ER = ER, TE = twin.E))
}


###========================================================================###
###============== LIST OF VERTEX COLORINGS AND EDGE COLORINGS =============###
###========================================================================###

### INPUT: 
#         p: a (even) numeric (number of vertices)
#         g: a list of (L.as, E, E.as)

### OUTPUT: a list of
#         v.list: a list of labels of vertex coloring 
#         e.list: a list of labels of edge coloring


rcox.lists <- function(p, g){
  L.as <- g$L.as
  E <- g$E
  E.as <- g$E.as
  q <- p/2
  
  v <- list()
  j <- 1
  for (i in 1:q) {
    if((i %in% L.as) == FALSE){
      v[[j]] <- list(vLabel(i, p), vLabel(tau(i, p), p))
      j <- j+1
    } else {
      v[[j]] <- list(vLabel(i, p))
      v[[j+1]] <- list(vLabel(tau(i, p), p))
      j <- j+2
    }
  }
  
  if(is.null(E)){
    e <- NULL
  } else {
    TE <- outEdges(E, p)$TE
    es <- matdiff(TE, E.as)$val # symmetric edges
    es.full <- rbind(es, tauMat(es, p))
    if(is.null(es)){
      e <- list()
    } else {
      es.lab <- apply(es.full, 2, vLabel, p)
      e <- lapply(seq_len(nrow(es)), function(i) {
        list(es.lab[i,], es.lab[i+nrow(es),] )
      })
    }
    
    ea.full <- matdiff(E, es.full)$val # asymmetric edges
    if(is.null(ea.full)){
      e <- e
    } else {
      if(nrow(ea.full) == 1) {
        ea.lab <- t(as.matrix(vLabel(ea.full[1,], p))) #if not do that, it return vector, not matrix
      } else {
        ea.lab <- apply(ea.full, 2, vLabel, p)
      }
      #l <- length(e)
      e <- append(e, lapply(seq_len(nrow(ea.full)), function(i) {list(ea.lab[i,])}))
    }
  }
  return(list(v.list = v, e.list = e))
}


###========================================================================###
###== P-VALUES FROM THE LIKELIHOOD RATIO TEST RELATIVE TO THE FULL MODEL ==###
###========================================================================###

### INPUT: 
#         p: a (even) numeric (number of vertices)
#         g: a list of (L.as, E, E.as)
#         logLik0: a value of log-likelihood of the full model obtained from gRc::rcox()
#         dof0: degree of freedom of the full model obtained from gRc::rcox()
#         S: an empirical covariance matrix (computed from data)
#         n: the number of observations

### OUTPUT: the p-value, which is the probability of obtaining a chi-square

### PACKAGE: gRc
# Hojsgaard, S., Lauritzen, S. L., et al. (2007). Inference in graphical Gaussian models 
# with edge and vertex symmetries with the gRc package for R. Journal of Statistical Software.


pval.function <- function(p, g, logLik0, dof0, S, n){
  VE.list <- rcox.lists(p, g)
  fit <- gRc::rcox(vcc = VE.list$v.list, ecc = VE.list$e.list,
                   S = S, n = n, method = "ipms")
  logLik <- fit$fitInfo$logL
  dof <- length(VE.list$e.list) + length(VE.list$v.list)
  
  teststat <- -2 * (logLik - logLik0)
  pval <- pchisq(teststat, df = dof0-dof, lower.tail = FALSE)
  return(pval)
}


###========================================================================###
###=========================== LAYOUT FUNCTION ============================###
###========================================================================###

### INPUT: 
#         p: a (even) numeric (number of vertices)

### OUTPUT: a two-column matrix where i-th row indicates the coordinates of i-th node,
#           for i = 1,...,p, and 1st column refer to x-coor, 2nd column refer to y-coor


layoutSym <- function(p, seed = 15){
  q <- p/2
  coor <- matrix(NA, ncol = 2, nrow = p)
  if(p == 4){
    coor[1:2, 1] <- coor[2,2] <- coor[4,2] <- 0
    coor[3:4, 1] <- coor[1,2] <- coor[3,2] <- 2
  } else {
    middle <- ifelse(p < 40, 6, 8)
    if(is.numeric(seed)) set.seed(seed)
    temp <- 1.2*runif(q,0, middle - 1.5)
    dis <- abs(temp-middle)
    for (i in 1:p) {
      if(i %in% 1:q){
        coor[i,1] <- temp[i]
        coor[i,2] <- -i/3
      } else {
        coor[i,1] <- temp[i-q] + 2*dis[i-q]
        coor[i,2] <- -(i-q)/3
      }
    }
  }
  return(coor)
}

###========================================================================###
###======================= VISUALIZATION FROM MODELS ======================###
###========================================================================###

### INPUT: 
#         p: a (even) numeric (number of vertices)
#         g: a list of (L.as, E, E.as)
#         seed: a single value, interpreted as an integer

### OUTPUT: plotting of colored graph

outGraph <- function(p, g, seed){
  L.as <- g$L.as
  E <- g$E
  E.as <- g$E.as
  q <- p/2
  
  # vertices
  L.color <- rep("white", q)
  temp <- setdiff(1:q, L.as)
  L.color[temp] <- rainbow(length(temp))
  V.color <- rep(L.color, 2) # color of vertices
  V.label <- vLabel(1:p, p) # label of vertices
  
  # edges
  if(!is.null(E)){
    
    # graphs
    adj.mat <- igraph::get.adjacency(graph.edgelist(E))
    adj.mat <- as.matrix(adj.mat)
    adj.mat <- adj.mat + t(adj.mat)
    g <- igraph::graph_from_adjacency_matrix(adj.mat, weighted = TRUE, 
                                             mode = "undirected", diag = FALSE)
    
    # integrate matrix E follow the edge list from g (i.e E(g))
    E.color <- rep("black", nrow(E))
    E <- igraph::get.edgelist(g, names = TRUE)
    
    twin.E <- outEdges(E, p)$TE
    Es <- matdiff(twin.E, E.as)$val # symmetric edges
    if(is.null(Es) == FALSE){
      Ecolors <- rainbow(nrow(Es))
      for (i in 1:nrow(Es)) {
        ind.colorL <- intersectMat(E, Es[i,, drop = FALSE])$ind.row
        ind.colorR <- intersectMat(E, tauMat(Es[i,, drop = FALSE],p))$ind.row
        E.color[ind.colorL] <- E.color[ind.colorR] <- Ecolors[i]
      }
    }
    plot.g <- plot(g, vertex.color = V.color, edge.color = E.color, vertex.label = V.label,
                   vertex.label.dist = 0, vertex.size = 15, layout = layoutSym(p, seed), edge.width = 3)
  } else {
    # graphs
    g <- igraph::make_empty_graph(n = p, directed = FALSE)
    plot.g <- plot(g, vertex.color = V.color, vertex.label = V.label, vertex.label.dist = 0,
                   vertex.size = 15, layout = layoutSym(p, seed))
  }
  return(plot.g)
}



