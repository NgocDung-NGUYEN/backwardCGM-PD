
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
