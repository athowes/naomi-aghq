# sparse_quad <- fit_aghq(tmb_inputs, k = 2, basegrid = sparse_grid, control = control)
# fit_aghq <- function(tmb_input, inner_verbose = FALSE, progress = NULL, map = NULL, ...)

tmb_input <- tmb_inputs
inner_verbose <- FALSE
progress <- NULL
map <- NULL

stopifnot(inherits(tmb_input, "naomi_tmb_input"))
obj <- local_make_tmb_obj(tmb_input$data, tmb_input$par_init, calc_outputs = 0L, inner_verbose, progress, map)

# Start of expose quad <- aghq::marginal_laplace_tmb(obj, startingvalue = obj$par, k = 2, basegrid = sparse_grid, control = control) ----

ff <- obj
k <- 2
startingvalue <- obj$par
transformation <- default_transformation()
optresults <- NULL
basegrid <- sparse_grid
control <- control

validate_control(control, type = "tmb")
validate_transformation(transformation)
transformation <- make_transformation(transformation)

thetanames <- NULL

if (exists("par", ff)) thetanames <- make.unique(names(ff$par), sep = "")

if (control$numhessian) {
  ff$he <- function(theta) numDeriv::jacobian(ff$gr, theta, method = "Richardson")
}

# Start of expose quad <- aghq(ff = ff, k = k, transformation = transformation, startingvalue = startingvalue, optresults = optresults, basegrid = basegrid, control = control) ----

validate_control(control)
validate_transformation(transformation)
transformation <- make_transformation(transformation)

if (is.null(optresults)) utils::capture.output(optresults <- optimize_theta(ff, startingvalue, control))

normalized_posterior <- normalize_logpost(optresults, k, basegrid = basegrid, ndConstruction = control$ndConstruction)

if (control$onlynormconst) return(normalized_posterior$lognormconst)

out <- list(normalized_posterior = normalized_posterior, optresults = optresults, control = control, transformation = transformation)

class(out) <- "aghq"
d <- length(startingvalue)
marginals <- vector(mode = "list", length = d)

# save.image(file = "env.rdata")

if (control$method_summaries[1] == "correct") {
  for (j in 1:d) marginals[[j]] <- aghq:::marginal_posterior.aghq(out, j, method = "correct")
} else {
  for (j in 1:d) marginals[[j]] <- aghq:::marginal_posterior.aghq(out, j, method = "reuse")
}

j <- 17
marginals[[j]] <- aghq:::marginal_posterior.aghq(out, j, method = "correct")

# Error in rescale.NIGrid(thegrid, m = m, C = Matrix::forceSymmetric(solve(H)),  :
#   (in rescale) no appropriate covariance matrix C

# Start of expose aghq:::marginal_posterior.aghq(out, j = 17, method = "correct") ----

quad <- out
j <- 17
qq <- NULL
method <- "correct"

method <- method[1]
if (method == "auto") method <- "reuse"
if (method == "reuse") return(marginal_posterior.list(quad$optresults, get_numquadpoints(quad), j))

S <- get_param_dim(quad)
idxorder <- c(j, (1:S)[-j])
thetaminusj <- (1:S)[-j]

mm <- get_mode(quad)[idxorder]
HH <- get_hessian(quad)[idxorder, idxorder]
gg <- mvQuad::createNIGrid(1, "GHe", get_numquadpoints(quad))
mvQuad::rescale(gg, m = mm[1], C = solve(HH)[1, 1], dec.type = 2)
qqq <- as.numeric(mvQuad::getNodes(gg))
out <- vector(mode = "list", length = length(qqq))

for (i in 1:length(qqq)) {
  out[[i]] <- aghq:::marginal_posterior.aghq(quad = quad, j = j, qq = qqq[i], method = "correct")
}

# Start of expose aghq:::marginal_posterior.aghq(quad = quad, j = j, qq = qqq[1], method = "correct") ----

i <- 1
quad <- quad
j <- j
qq <- qqq[i]
method <- "correct"

method <- method[1]

S <- get_param_dim(quad)
idxorder <- c(j, (1:S)[-j])
thetaminusj <- (1:S)[-j]

cname <- colnames(get_nodesandweights(quad))[j]

fn <- function(theta) quad$optresults$ff$fn(aghq:::splice(theta, qq, j))
gr <- function(theta) quad$optresults$ff$gr(aghq:::splice(theta, qq, j))[-j]
he <- function(theta) quad$optresults$ff$he(aghq:::splice(theta, qq, j))[-j, -j]

ffm <- list(fn = fn, gr = gr, he = he)
newcontrol <- quad$control
newcontrol$onlynormconst <- TRUE
newcontrol$negate <- FALSE

lognumerator <- get_log_normconst(aghq(ffm, get_numquadpoints(quad), get_mode(quad)[-j], control = newcontrol))

# Start expose quadm <- aghq(ffm, k = aghq::get_numquadpoints(quad), aghq::get_mode(quad)[-j], control = newcontrol)

ff <- ffm
k <- aghq::get_numquadpoints(quad)
startingvalue <- aghq::get_mode(quad)[-j]
control <- newcontrol
basegrid <- NULL
optresults <- NULL
transformation <- default_transformation()

validate_control(control)
validate_transformation(transformation)
transformation <- make_transformation(transformation)

if (is.null(optresults)) utils::capture.output(optresults <- aghq::optimize_theta(ff, startingvalue, control))

normalized_posterior <- normalize_logpost(optresults, k, basegrid = basegrid, ndConstruction = control$ndConstruction)

# Start expose normalized_posterior <- normalize_logpost(optresults, k, basegrid = basegrid, ndConstruction = control$ndConstruction)

optresults <- optresults
whichfirst <- 1
k <- k
basegrid <- basegrid
ndConstruction <- control$ndConstruction

S <- length(optresults$mode)

thegrid <- mvQuad::createNIGrid(dim = S, type = "GHe", level = k, ndConstruction = ndConstruction)

idxorder <- c(whichfirst, (1:S)[-whichfirst])

m <- optresults$mode[idxorder]
H <- optresults$hessian[idxorder, idxorder]

mvQuad::rescale(thegrid, m = m, C = Matrix::forceSymmetric(solve(H)), dec.type = 2)

#' Error is here!
#' H is not symmetric, as evidenced by imaginary eigenvalues
eigen(H, only.values = TRUE)

#' The asymmetry is relatively isolated
skew_part <- H - Matrix::forceSymmetric(H)
image(skew_part)

sum(skew_part)
plot(abs(skew_part)[order(abs(skew_part))])

#' Which parameter name does theta17 correspond to?
thetanames[17]

#' Which elements of the Hessian are very asymmetric? Look at the 10 largest
asymmetric_indices <- head(order(abs(skew_part), decreasing = TRUE), n = 10)
skew_part[asymmetric_indices]

data.frame(row = asymmetric_indices %% nrow(H), col = asymmetric_indices %/% nrow(H) + 1)

#' What is the 10th parameter where the greatest asymmetries are?
thetanames[10]

if (control$onlynormconst) quadm <- normalized_posterior$lognormconst

lognumerator <- get_log_normconst(quadm)
out <- data.frame(qq, lognumerator - get_log_normconst(quad))
colnames(out) <- c(cname, "logmargpost")

# End of expose aghq:::marginal_posterior.aghq(quad = quad, j = j, qq = qqq[1], method = "correct") ----

out <- Reduce(rbind, out)

out

# End of expose quad <- aghq(ff = ff, k = k, transformation = transformation, startingvalue = startingvalue, optresults = optresults, basegrid = basegrid, control = control) ----

if (control$onlynormconst) return(quad)

distinctthetas <- quad$normalized_posterior$nodesandweights[, grep("theta", colnames(quad$normalized_posterior$nodesandweights))]

if (!is.data.frame(distinctthetas))
  distinctthetas <- data.frame(theta1 = distinctthetas)
modesandhessians <- distinctthetas
if (is.null(thetanames)) {
  thetanames <- colnames(distinctthetas)
} else {
  colnames(modesandhessians) <- thetanames
  colnames(quad$normalized_posterior$nodesandweights)[grep("theta", colnames(quad$normalized_posterior$nodesandweights))] <- thetanames
}

modesandhessians$mode <- vector(mode = "list", length = nrow(distinctthetas))
modesandhessians$H <- vector(mode = "list", length = nrow(distinctthetas))

for (i in 1:nrow(distinctthetas)) {
  theta <- as.numeric(modesandhessians[i, thetanames])
  ff$fn(theta)
  mm <- ff$env$last.par
  modesandhessians[i, "mode"] <- list(list(mm[ff$env$random]))
  H <- ff$env$spHess(mm, random = TRUE)
  H <- rlang::duplicate(H)
  modesandhessians[i, "H"] <- list(list(H))
}

quad$modesandhessians <- modesandhessians
class(quad) <- c("marginallaplace", "aghq")

quad$obj <- obj
quad

# End of expose quad <- aghq::marginal_laplace_tmb(obj, startingvalue = obj$par, k = 2, basegrid = sparse_grid, control = control) ----
