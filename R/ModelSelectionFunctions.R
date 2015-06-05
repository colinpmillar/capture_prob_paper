phiBIC <- function(object, ..., phi = 2.106) {
  ll <- if ("stats4" %in% loadedNamespaces()) stats4:::logLik
    else logLik
  Nobs <- if ("stats4" %in% loadedNamespaces()) stats4:::nobs
    else nobs
  lls <- ll(object)
  nos <- attr(lls, "nobs")
  -2 * as.numeric(lls)/phi + log(nos) * attr(lls, "df")
}

# should maybe return the best model as well summaries...
summaryMods <- function(lst, m0 = NULL, order = TRUE, fn = phiBIC, phi = 2.11) {
  #aics <- sapply(lst, AIC)
  fn <- match.fun(fn)
  aics <- sapply(lst, fn, phi = phi)
  lliks <- sapply(lst, logLik) / phi
  dfs <- sapply(lst, function(x) attr(logLik(x), "df"))

  tab <- 
   data.frame(
    forms = sapply(lst, function(x) paste(deparse(x$formula, width.cutoff = 500L))),
    aic = aics,
    llik = lliks,
    df = dfs,
    phi = phi
    )

  if (!is.null(m0)) {
    tab $ Daic <- tab $ aic - fn(m0, phi = phi)
    tab $ Chisqp <- 1-pchisq(abs(lliks - logLik(m0)/phi), 
                           abs(attr(logLik(m0), "df") - dfs))
  }
  if (order) tab <- tab[order(aics),]

  unique(tab)  
}

runModels <- function(chosen, f1s, data = ef3, fn = phiBIC, phi = 2.11, ...) {
  # build a formula list of additions
  formsadd <- lapply(f1s[!chosen], function(x) as.formula(paste0("X ~ ", paste(c(f1s[chosen], x), collapse = " + "))))

  descadd <- f1s[!chosen]
  dropadd <- rep("add", sum(!chosen))

  #build a formula list dropping elements
  if (any(chosen)) {
    if (sum(chosen)==1) {
      formsdrop <- list(X ~ 1)
    } else {
      formsdrop <- lapply(seq(sum(chosen)), function(i) as.formula(paste0("X ~ ", paste(f1s[chosen][-i], collapse = " + ") )))
    }
    descdrop <- f1s[chosen]
    forms <- c(formsadd, formsdrop)
    desc <- c(descadd, descdrop)
    dropadd <- c(dropadd, rep("drop", sum(chosen)))
  } else {
    forms <- formsadd
    desc <- descadd
  }

  mods <- lapply(forms, efp, data = data, passes = "Runs", ...)

  if (all(!chosen)) {
    m0 <- efp(X ~ 1, data = data, passes = "Runs", ...)
  } else {
    m0 <- efp(as.formula(paste0("X ~ ", paste(f1s[chosen], collapse = " + "))), data = data, passes = "Runs", ...)    
  }
  
  tab <- cbind(what = desc, step = dropadd, summaryMods(mods, m0 = m0, order = FALSE, fn = fn, phi = phi)[,-1] )
  tab[order(tab $ Daic),]
}
  
runSelection <- function(forms, data = ef3, fn = phiBIC, phi = 2.11, start = NULL) {
  chosen <- rep(FALSE, length(forms))
  if (!is.null(start)) chosen <- start
  out <- list(Daic = -1)
  tol <- 0
  outlist <- list() # grow this - bad!
  i <- 0
  while(out $ Daic[1] < tol & !all(chosen)) {
    out <- runModels(chosen, forms, verbose = FALSE, data = data, fn = fn, phi = phi)
    print(head(out, 10), digits = 3)
    outlist[[i <- i + 1]] <- out
    if ( out $ Daic[1] < 0) {
      # then model is an improvement
      which <- which(forms %in% out $ what[1])
      chosen[which] <- !chosen[which]
      cat("\t", paste(forms[chosen], collapse = " + "), "\n ----------------------------------------------- \n")
    } 
  }

  if (all(!chosen)) {
    mod <- efp(X ~ 1, data = data, passes = "Runs")
  } else {
    mod <- efp(as.formula(paste0("X ~ ", paste(forms[chosen], collapse = " + "))), data = data, passes = "Runs")
  }

  list(history = outlist, chosen = chosen, mod = mod)
}


getScale <- function(wk, mod) {
  # get p predictions for ef
  wk $ pbig <- fitted(mod)

  # get logLik components
  wk $ R <- with(wk, s - 1 - Z)
  wk $ llsat <- with(wk, T * log(psat) + T * R * log(1-psat) - T * log(1 - (1-psat)^s) )
  wk $ llbig <- with(wk, T * log(pbig) + T * R * log(1-pbig) - T * log(1 - (1-pbig)^s) )

  # calculate Deviance components and residuals
  wk $ devcomp <- with(wk, 2 * (llsat - llbig))
  wk $ devres <- with(wk, sign(psat - pbig)*sqrt(devcomp))

  # scale estimate
  sum(wk $ devcomp)/(nrow(wk)-mod $ rank)
}
