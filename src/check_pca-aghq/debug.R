createNIGrid <- function(dim=NULL, type=NULL, level=NULL, ndConstruction="product", level.trans=NULL){

  if (is.null(dim) | (dim%%1)!=0 | dim < 1) {
    stop("'dim' is not appropriate defined")
  }

  if (class(type)!="list"){
    if (class(type)=="character") type <- as.list(type)
    if (class(type)=="costumRule") type <- list(type)
  }

  if (is.null(type) | (length(type) > 1 & length(type) < dim)) {
    stop("'type' is not appropriate defined")
  }

  if (length(type) < dim){
    type <- replicate(type[[1]], n=dim, FALSE)
  }

  if (is.null(level) | any(level%%1!=0)){
    stop("'level' is not appropriate defined")
  }

  if (length(level) < dim){
    level <- rep(level, dim)
  }

  level <- matrix(level, ncol = dim)

  if (is.logical(level.trans)) {
    if (level.trans == TRUE){
      level.trans <- function(x){2^(x-1)}
    }
  }

  if (!is.function(level.trans)){
    level.trans <- function(x){x}
  }

  if (dim > 1){
    if (ndConstruction=="product"){
      nw <- .fgridnD(dim = dim, type = type, level = level.trans(level))
    }
    if (ndConstruction=="sparse"){
      nw <- .sgridnD(dim = dim, type = type, level = level[1], level.trans = level.trans)
    }
  } else {
    nw <- .grid1D(type = type[[1]], level = level.trans(level))
  }

  object=new.env(parent=globalenv())

  object$dim <- dim

  type.text <- lapply(type,
                      function(tmp){
                        if (class(tmp)=="costumRule") return("costum")
                        if (class(tmp)=="character")  return(tmp)
                      })
  object$type <- unlist(type.text)
  object$level <- level
  object$level.trans <- level.trans

  object$ndConstruction <- ndConstruction

  if (exists("features", where=nw)){
    object$features <- c(type="static", move=0L, nw[["features"]])
  } else {
    object$features <- list(type="static", move=0L, initial.domain = matrix(c(NA, NA), ncol=2))
  }

  if (!exists("n", where=nw) | !exists("n", where=nw)) stop("this quadrature rule doesn't provide nodes (or/and) weights")

  object$nodes <- as.matrix(nw[["n"]])
  object$weights <- as.matrix(nw[["w"]])
  class(object) <- "NIGrid"
  return(object)
}


.fgridnD <- function(dim, type, level) {
  le <- numeric(dim)
  domain <- matrix(NA, nrow = dim, ncol = 2)
  storage <- vector(mode = "list", length = dim)
  for (i in 1:dim) {
    tmp <- .grid1D(type[[i]], level[i])
    le[i] <- length(tmp$w)
    domain[i, ] <- tmp$features$initial.domain
    storage[[i]] <- tmp
  }

  No.Points <- prod(le)
  w <- rep(1, No.Points)
  n <- matrix(NA, nrow = No.Points, ncol = dim)

  for (i in 1 :dim) {
    storage[[i]]
    n[ ,i] = rep(storage[[i]]$n, each = prod(le[1:i - 1]), times = (No.Points/prod(le[1:i])))
    w = w * rep(storage[[i]]$w, each = prod(le[1:i - 1]), times = (No.Points/prod(le[1:i])))
  }
  return(list(n = n, w = w, features = list(initial.domain = domain)))
}

.hardCodedTypes <- c("cNC1", "cNC2", "cNC3", "cNC4", "cNC5", "cNC6",
                     "oNC0", "oNC1", "oNC2", "oNC3")

.extGaussQuad <- c("GLe", "GLa", "GHe", "GHN")

.preDefinedTypes <- structure(list(nLe = structure(list(levels = c(1L, 3L, 3L, 7L, 7L, 7L, 15L,
                                                                   15L, 15L, 15L, 15L, 15L, 31L, 31L, 31L, 31L, 31L, 31L, 31L,
                                                                   31L, 31L, 31L, 31L, 31L, 63L)), .Names = "levels"),
                                   GKr = structure(list(levels = c(NA, NA, 3L, NA, 5L, NA, 7L, NA, 9L, NA, 11L,
                                                                   NA, 13L, NA, 15L, NA, 17L, NA, 19L, NA, 21L, NA, 23L,
                                                                   NA, 25L, NA, 27L, NA, 29L)), .Names = "levels"),
                                   nHe = structure(list(levels = c(1L, 3L, 3L, 7L, 9L, 9L, 9L, 9L, 17L, 19L,
                                                                   19L, 19L, 19L, 19L, 19L, 31L, 33L, 35L, 35L, 35L, 35L,
                                                                   35L, 35L, 35L, 35L)), .Names = "levels"),
                                   nHN = structure(list(levels = c(1L, 3L, 3L, 7L, 9L, 9L, 9L, 9L, 17L, 19L,
                                                                   19L, 19L, 19L, 19L, 19L, 31L, 33L, 35L, 35L, 35L, 35L,
                                                                   35L, 35L, 35L, 35L)), .Names = "levels"),
                                   Leja = structure(list(levels = 1:141), .Names = "levels")),
                              .Names = c("nLe", "GKr", "nHe", "nHN", "Leja"))

.grid1D <- function(type, level=1) {

  if (class(type)=="costumRule"){
    if (!(level %in% c(1:length(type)))){
      stop(paste("degree (",level,") for user defined rule not supported \n") )
    } else {
      return(type[[level]])
    }
  }

  if (class(type)!="character") stop("type of quadrature rule not appropriate defined")

  if (!((type %in% c(.hardCodedTypes, names(.preDefinedTypes), .extGaussQuad)))) {
    if (!existsFunction(type)) {
      stop("type of quadrature rule is not supported")
    } else {
      tmp.fun <- match.fun(type)
      return(tmp.fun(level))
    }
  }

  if (type %in% .hardCodedTypes) {
    quad.type <- substr(type,1,3)

    # closed-Newton-Cotes forumla
    if (quad.type=="cNC") {

      degree <- as.numeric(substr(type,4,4))
      if (!(degree %in% c(1:6)) ) stop("degree of closed-Newton-Cotes formula must between 1 and 6")

      if (degree == 1) weights <- c(1/2, 1/2)
      if (degree == 2) weights <- c(1/6, 4/6, 1/6)
      if (degree == 3) weights <- c(1/8, 3/8, 3/8, 1/8)
      if (degree == 4) weights <- c(7/90, 32/90, 12/90, 32/90, 7/90)
      if (degree == 5) weights <- c(19/288, 75/288, 50/288, 50/288, 75/288, 19/288)
      if (degree == 6) weights <- c(41/840, 216/840, 27/840, 272/840, 27/840, 216/840, 41/840)

      interv_length <- 1/level
      weights <- weights * interv_length
      steps <- interv_length/degree
      n <- seq(0,1,steps)
      w <- numeric(length(n))

      for (i in 1:level) {
        w[((i-1)*degree+1):(i*degree+1)] <- w[((i-1)*degree+1):(i*degree+1)] + weights
      }

      return(list(n = n,w = w, features=list(initial.domain=matrix(c(0,1), ncol = 2))))
    }

    # open-Newton-Cotes Forumla
    if (quad.type=="oNC") {
      degree <- as.numeric(substr(type,4,4))
      if (!(degree %in% c(0:3)) ) stop("degree of open-Newton-Cotes formula must between 0 and 3")

      if (degree == 0) {
        weights <- c(1)
        nodes <- c(1/2)
      }
      if (degree == 1) {
        weights <- c(1/2, 1/2)
        nodes <- c(1/4, 3/4)
      }
      if (degree == 2) {
        weights <- c(3/8, 2/8, 3/8)
        nodes <- c(1/6, 1/2, 5/6)
      }
      if (degree == 3) {
        weights <- c(13/48, 11/48, 11/48, 13/48)
        nodes <- c(1/8, 3/8, 5/8, 7/8)
      }
      interv_length <- 1/level
      weights <- weights * interv_length
      nodes <- nodes * interv_length

      n = numeric((degree+1)*level)
      w = rep(weights,level)

      for (i in 1:level){
        n[((i-1)*(degree+1)+1):(i*(degree+1))] <- (i-1)*interv_length + nodes
      }

      return(list(n = n, w = w, features=list(initial.domain=matrix(c(0,1), ncol = 2))))
    }
  }

  if (type %in% names(.preDefinedTypes)) {
    if (!(level %in% c(1:length(QuadRules[[type]]))) ){
      stop(paste("degree (",level,") for rule '", type, "' not supported \n") )
    } else {
      return(QuadRules[[type]][[level]])
    }
  }

  if (type %in% .extGaussQuad) {
    if (type=="GLe") {
      tmp.nw <- statmod::gauss.quad(level, kind="legendre")
      tmp.nw$nodes <- tmp.nw$nodes / 2 + 0.5
      tmp.nw$weights <- tmp.nw$weights / 2
      tmp.dom <- list(initial.domain=matrix(c(0,1), ncol = 2))
    }

    if (type=="GLa") {
      tmp.nw <- statmod::gauss.quad(level, kind="laguerre")
      tmp.dom <- list(initial.domain=matrix(c(0,Inf), ncol = 2))
    }

    if (type=="GHe") {
      tmp.nw <- statmod::gauss.quad(level, kind="hermite")
      tmp.nw$nodes <- tmp.nw$nodes * sqrt(2)
      tmp.nw$weights <- tmp.nw$weights * sqrt(2) / sqrt(exp(-(tmp.nw$nodes^2)))
      tmp.dom <- list(initial.domain=matrix(c(-Inf, Inf), ncol = 2))
    }
    if (type=="GHN") {
      tmp.nw <- statmod::gauss.quad(level, kind="hermite")
      tmp.nw$nodes <- tmp.nw$nodes * sqrt(2)
      tmp.nw$weights <- tmp.nw$weights * sqrt(2) / sqrt(exp(-(tmp.nw$nodes^2)))*dnorm(tmp.nw$nodes)
      tmp.dom <- list(initial.domain=matrix(c(-Inf, Inf), ncol = 2))
    }

    return(list(n = as.numeric(tmp.nw$nodes), w = as.numeric(tmp.nw$weights), features=tmp.dom))
  }

}

.grid1D(type = "GHe", level = 3)
result <- createNIGrid(dim = 2, type = "GHe", level = c(3, 1))

#' Expose createNIGrid
dim <- 2
type <- "GHe"
level <- c(3, 1)
level.trans <- NULL

if (length(type) < dim){
  type <- replicate(type[[1]], n=dim, FALSE)
}

if (length(level) < dim){
  level <- rep(level, dim)
}

level <- matrix(level, ncol = dim)

if (is.logical(level.trans)) {
  if (level.trans == TRUE){
    level.trans <- function(x){2^(x-1)}
  }
}

if (!is.function(level.trans)){
  level.trans <- function(x){x}
}

#' This is the place the weights get calculated!
nw <- .fgridnD(dim = dim, type = type, level = level.trans(level))

#' Expose .fgridnD
dim <- dim
type <- type
level <- level.trans(level)

le <- numeric(dim)
domain <- matrix(NA, nrow = dim, ncol = 2)
storage <- vector(mode = "list", length = dim)
for (i in 1:dim) {
  tmp <- .grid1D(type[[i]], level[i])
  le[i] <- length(tmp$w)
  domain[i, ] <- tmp$features$initial.domain
  storage[[i]] <- tmp
}

No.Points <- prod(le)
w <- rep(1, No.Points)
n <- matrix(NA, nrow = No.Points, ncol = dim)

#' Break a loop over i = 1, 2 down
storage[[1]]
n[ ,1] = rep(storage[[1]]$n, each = prod(le[1:1 - 1]), times = (No.Points/prod(le[1:1])))
w = w * rep(storage[[1]]$w, each = prod(le[1:1 - 1]), times = (No.Points/prod(le[1:1])))

storage[[2]]
n[ ,2] = rep(storage[[2]]$n, each = prod(le[1:2 - 1]), times = (No.Points/prod(le[1:2])))
w = w * rep(storage[[2]]$w, each = prod(le[1:2 - 1]), times = (No.Points/prod(le[1:2])))

list(n = n, w = w, features = list(initial.domain = domain))

mvQuad::getWeights(result)
