#'
#'  simulateCSCP.R
#'

require(spatstat)
require(spatstat.utils)

simulateCSCP <- function(p, W = owin(), nsim=1, drop=TRUE) {
  if(!is.list(p) || !all(sapply(p, is.list)))
    stop("Argument p should be a list of lists of parameter values",
         call.=FALSE)
  W <- as.owin(W)
  check.1.integer(nsim)
  drop <- isTRUE(drop)
  argh <- list(W=W, nsim=nsim, drop=FALSE) # sic
  ## 
  k <- length(p)
  Zlist <- vector(mode="list", length=k)
  for(i in seq_len(k))
    Zlist[[i]] <- do.call(simulateGRF, append(p[[i]], argh))
  Lamlist <- vector(mode="list", length=nsim)
  Ylist <- vector(mode="list", length=nsim)
  for(j in seq_len(nsim)) {
    Lamj <- 0
    for(i in seq_len(k)) {
      Zij <- Zlist[[i]][[j]]
      Lamj <- Lamj + Zij^2
    }
    Lamlist[[j]] <- Lamj
    Ylist[[j]] <- rpoispp(Lamj)
  }
  result <- simulationresult(Ylist, nsim, drop)
  attr(result, "Lambda") <- simulationresult(Lamlist, nsim, drop)
  return(result)
}

simulateGRF <- function(model=c("exponential", "gauss", "stable",
                                "gencauchy", "matern"),
                        ..., 
                        W=owin(), nsim=1, drop=FALSE) {
  model <- match.arg(model)
  switch(model,
         exponential = {
           rGRFexpo(W=W, ..., nsim=nsim, drop=drop)
         }, gauss = {
           rGRFgauss(W=W, ..., nsim=nsim, drop=drop)
         }, stable = {
           rGRFstable(W=W, ..., nsim=nsim, drop=drop)
         }, gencauchy = {
           rGRFgencauchy(W=W, ..., nsim=nsim, drop=drop)
         }, matern = {
           rGRFmatern(W=W, ..., nsim=nsim, drop=drop)
         }, stop(paste("Model", sQuote(model), "not recognised"),
                 call.=FALSE))
}




## --- Reparameterised wrapper: stationary CSCP with (barlambda, w, alphas)
## Minimal-change design:
##  - Convert (barlambda, w, alphas, scales) -> list-of-lists structure: list(list(mu=..., scale=..., var=...), ...)
##  - Calls existing simulateCSCP()
##  - If lambda0 > 0, add homogeneous Poisson background and superimpose

## Reparameterised wrapper (stationary): uses (barlambda, w, alphas, scales)
## Minimal change: no 'win' passed into simulateCSCP(); 'win' only used
## to generate the baseline Poisson and to set/confirm the window.

simulateCSCP_reparam <- function(barlambda,
                                 w,
                                 alphas,
                                 scales,
                                 nsim = 1,
                                 win  = NULL,
                                 ...) {
  stopifnot(is.numeric(barlambda), length(barlambda) == 1, barlambda > 0)
  stopifnot(is.numeric(w), length(w) == 1, w >= 0, w <= 1)
  stopifnot(is.numeric(alphas), all(alphas >= 0), sum(alphas) > 0)
  stopifnot(is.numeric(scales), length(scales) == length(alphas), all(scales > 0))
  
  ## normalise alphas
  alphas <- as.numeric(alphas) / sum(alphas)
  
  ## stationary mapping: sigma_i^2 = w * barlambda * alpha_i ; mu_i = 0
  vars  <- w * barlambda * alphas
  
  ## structure expected by simulateCSCP(): list of lists
  comps <- lapply(seq_along(alphas), function(i) {
    list(mu = 0, scale = scales[i], var = vars[i])
  })
  
  ## ---- call existing simulator (NO 'win' argument)
  cs_out <- simulateCSCP(comps, nsim = nsim, ...)
  
  ## baseline Poisson with intensity lambda0 = (1 - w) * barlambda
  lambda0 <- (1 - w) * barlambda
  if (lambda0 <= 0) return(cs_out)
  
  ## Determine window for adding the baseline component
  getW <- function(obj) {
    if (!is.null(win)) return(win)
    if (inherits(obj, "ppp")) return(Window(obj))
    if (inherits(obj, "solist") && length(obj) >= 1 && inherits(obj[[1]], "ppp"))
      return(Window(obj[[1]]))
    if (is.list(obj) && length(obj) >= 1 && inherits(obj[[1]], "ppp"))
      return(Window(obj[[1]]))
    stop("Unable to determine window; pass 'win=' explicitly.")
  }
  W <- getW(cs_out)
  
  ## superimpose Poisson background with intensity lambda0
  add_bg <- function(pp) {
    bg <- rpoispp(lambda0, win = W)
    superimpose(pp, bg, W = W)
  }
  
  if (inherits(cs_out, "ppp")) {
    return(add_bg(cs_out))
  } else if (inherits(cs_out, "solist")) {
    ## keep 'solist' class so plot() behaves like your original
    out_list <- lapply(cs_out, add_bg)
    return(do.call(solist, out_list))
  } else if (is.list(cs_out) && length(cs_out) >= 1 && inherits(cs_out[[1]], "ppp")) {
    ## attempt to preserve solist-like behaviour
    out_list <- lapply(cs_out, add_bg)
    return(do.call(solist, out_list))
  } else {
    stop("simulateCSCP() returned an unsupported type; expected ppp or (so)list of ppp.")
  }
}