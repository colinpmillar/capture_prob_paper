moran.negbin <- function (link = "logit") {
  # set up link function
  linktemp <- substitute(link)
  if (!is.character(linktemp)) linktemp <- deparse(linktemp)
  okLinks <- c("logit", "probit", "cloglog", "cauchit", "log")
  if (linktemp %in% okLinks) {
      stats <- make.link(linktemp)
  } else if (is.character(link)) {
      stats <- make.link(link)
      linktemp <- link
  } else {
    if (inherits(link, "link-glm")) {
      stats <- link
      if (!is.null(stats$name)) linktemp <- stats$name
    } else {
      stop(gettextf("link \"%s\" not available for moran.negbin family; available links are %s", 
              linktemp, paste(sQuote(okLinks), collapse = ", ")), domain = NA)
      }
  }

  # define functions for variance, mean deviance and AIC
  variance <- function(mu) mu * (1 - mu)

  validmu <- function(mu) all(is.finite(mu)) && all(mu > 0 & mu < 1)
  # or maybe:
  # validmu <- function(mu) all(mu > 0)
  valideta <- function(eta) TRUE
  #dev.resids <- function(y, mu, wt) .Call(C_binomial_dev_resids, y, mu, wt)

  #aic <- function(y, n, mu, wt, dev) {
  #  m <- if (any(n > 1)) n else wt
  #  -2 * sum(ifelse(m > 0, (wt/m), 0) * dbinom(round(m * y), round(m), mu, log = TRUE))
  #}


  initialize <- expression({
    if (any(y < 0)) stop("negative values not allowed for the moran negative binomial family")
    n <- rep(1, nobs)
    mustart <- y + (y == 0)/6
  })


  simfun <- function(object, nsim) {
    ftd <- fitted(object)
    warning("simulate function not implemented yet.")
    NA
    #rnegbin(nsim * length(ftd), ftd, .Theta)
  }

  structure(list(family = "moran.negbin", link = linktemp, linkfun = stats$linkfun, 
        linkinv = stats$linkinv, variance = variance, dev.resids = dev.resids, 
        aic = aic, mu.eta = stats$mu.eta, initialize = initialize, 
        validmu = validmu, valideta = stats$valideta, simulate = simfun), 
        class = "family")
}


    variance <- function(mu) mu + mu^2/.Theta

    dev.resids <- function(y, mu, wt) 2 * wt * (y * log(pmax(1, 
        y)/mu) - (y + .Theta) * log((y + .Theta)/(mu + .Theta)))

    aic <- function(y, n, mu, wt, dev) {
        term <- (y + .Theta) * log(mu + .Theta) - y * log(mu) + 
            lgamma(y + 1) - .Theta * log(.Theta) + lgamma(.Theta) - 
            lgamma(.Theta + y)
        2 * sum(term * wt)
    }






