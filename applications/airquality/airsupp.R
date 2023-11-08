
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

###========================================================================###
###================= CONVERT FROM (|L, E, |E) TO MATRIX ===================###
###========================================================================###

formation.matrix <- function(p, pdsets, pdmat){
  q <- p/2
  if(is.null(pdmat)){
    pdmat <- diag(p)
    L.as <- pdsets$L.as
    E <- pdsets$E
    E.as <- pdsets$E.as
    
    # for diagonal entries
    colored.vertex <- setdiff(1:q, L.as)
    if(length(colored.vertex) > 0){
      diag(pdmat)[c(colored.vertex, tau(colored.vertex, p))] <- 2
    }
    
    
    # for off-diagonal entries
    if(!is.null(E)){
      pdmat[E] <- 1
      TE <- outEdges(E, p)$TE
      Ea <- matdiff(TE, E.as)$val
      colored.edge <- rbind(Ea, tauMat(Ea, p))
      pdmat[colored.edge] <- 2 
      pdmat <-  pdmat + t(pdmat) - diag(diag(pdmat))
    }
    return(pdmat)
  }
  
  if(is.null(pdsets)){
    pdmat.LL <- pdmat[1:q, 1:q]
    L.as <- which(diag(pdmat.LL) == 1)
    if(length(L.as) == 0) {L.as <- NULL}
    
    pdmat[lower.tri(pdmat, diag = TRUE)] <- 0
    E <- which(pdmat != 0, arr.ind = T)
    if (length(E) == 0){
      E <- NULL
    } else { colnames(E) <- NULL }
    
    
    full.Eas <- which(pdmat == 1, arr.ind = T)
    colnames(full.Eas) <- NULL
    E.as <- outEdges(full.Eas, p)$TE
    
    pdsets <- list(L.as = L.as, E = E, E.as = E.as)
    return(pdsets)
  }
}


###=============================================================###
###================= PLOT PD-CGs from MATRIX ===================###
###=============================================================###

plotG <- function(matrix, p, pvalue){
  twos <- which(matrix == 2, arr.ind = TRUE)
  twos <- twos[twos[,1] >= twos[,2],]
  
  zeros <- which(matrix == 0, arr.ind = TRUE)
  zeros <- zeros[zeros[,1] > zeros[,2],]
  
  ones <- which(matrix == 1, arr.ind = TRUE)
  ones.sm <- ones[ones[,1] < ones[,2],]
  ones.ve <- ones[ones[,1] == ones[,2],]
  
  colnames(ones.ve) <- colnames(ones.sm) <- colnames(zeros) <- colnames(twos) <- NULL
  ones.twin <- outEdges(ones.sm, dim(data.res)[2])$TE
  ones.twin <- rbind(ones.twin, tauMat(ones.twin, p))
  ones.twin <- rbind(ones.twin, ones.ve)
  
  only.ones <- matdiff(ones.sm, ones.twin)$val
  
  only.ones[ , c(1,2)] <- only.ones[ , c(2,1)] 
  ones.twin[ , c(1,2)] <- ones.twin[ , c(2,1)] 
  
  if(!is.null(only.ones) | length(only.ones) > 0){
    only.ones[,2] <- -only.ones[,2]
  }
  
  if(!is.null(ones.twin) | length(ones.twin) > 0){
    ones.twin[,2] <- -ones.twin[,2]
  }
  
  if(!is.null(twos) | length(twos) > 0){
    twos[,2] <- -twos[,2]
  }
  
  #Plotting
  par(mar=c(0,0,1,0), bty="n", pch="O", lty=2)
  
  g <- plot(c(-1,ncol(matrix)), c(0,-nrow(matrix)-1), type="n", 
            xlab=NA, ylab=NA, asp=1)
  
  lines(c((ncol(matrix)+1)/2,(ncol(matrix)+1)/2), 
        c(nrow(matrix),-nrow(matrix)-1), 
        type = "l", lty=1, col = "black")
  
  lines(c(0, ncol(matrix)+1), 
        c(-(nrow(matrix)+1)/2, -(nrow(matrix)+1)/2), 
        type = "l", lty=1, col = "black")
  
  points(only.ones, pch = 15, cex = 3, col = "black")
  points(ones.twin, pch = 19, cex = 3, col = "black")
  points(twos, pch = 21, cex = 3, col = "black", bg = "lightgray")
  
  if(!is.null(pvalue)){
    text(2, -nrow(matrix) + 2, paste0("pvalue = ", pvalue))
  }
  
  return(g)
}

###============================================================================###
###======== COMPUTE pvalue WITH INITIAL VALUES OF CONCENTRATION MATRIX ========###
###============================================================================###

pval.function <- function(p, g, logLik0, dof0, S, n, Kstart){
  VE.list <- rcox.lists(p, g)
  fit <- gRc::rcox(vcc = VE.list$v.list, ecc = VE.list$e.list, Kstart = Kstart,
                   S = S, n = n, method = "ipms")
  logLik <- fit$fitInfo$logL
  dof <- length(VE.list$e.list) + length(VE.list$v.list)
  
  teststat <- -2 * (logLik - logLik0)
  pval <- pchisq(teststat, df = dof0-dof, lower.tail = FALSE)
  return(pval)
}
