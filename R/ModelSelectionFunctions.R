phiBIC <- function(object, ..., phi = 2.106, nobs = NULL) {
  ll <- if ("stats4" %in% loadedNamespaces()) stats4:::logLik
    else logLik
  lls <- ll(object)
  nos <- if (is.null(nobs)) attr(lls, "nobs") else nobs
  -2 * as.numeric(lls)/phi + log(nos) * attr(lls, "df")
}

# should maybe return the best model as well summaries...
summaryMods <- function(lst, m0 = NULL, order = TRUE, fn = phiBIC, phi = 2.891, nobs = NULL) {
  #aics <- sapply(lst, AIC)
  fn <- match.fun(fn)
  ics <- sapply(lst, fn, phi = phi, nobs = nobs)
  lliks <- sapply(lst, logLik) / phi
  dfs <- sapply(lst, function(x) attr(logLik(x), "df"))
  Nobs <- sapply(lst, function(x) attr(logLik(x), "nobs"))
  if (!is.null(nobs)) Nobs[] <- nobs

  tab <- 
   data.frame(
    forms = sapply(lst, function(x) paste(deparse(x$formula, width.cutoff = 500L))),
    ic = ics,
    llik = lliks,
    df = dfs,
    phi = phi
    )

  if (!is.null(m0)) {
    tab $ Dic <- tab $ ic - fn(m0, phi = phi, nobs = nobs)

    kdiff <- abs(attr(logLik(m0), "df") - dfs)

    tab $ Fp <- 1 - pf(2*abs(lliks - logLik(m0)/phi)/kdiff, 
                           kdiff,
                           nobs - pmax(attr(logLik(m0), "df"), dfs) - 1)

  }
  if (order) tab <- tab[order(ics),]

  unique(tab)  
}

runModels <- function(chosen, forms1, data = ef, fn = phiBIC, phi = 2.981, nobs) {
  # build a formula list of additions
  formsadd <- lapply(forms1[!chosen], function(x) as.formula(paste0("n ~ ", paste(c(forms1[chosen], x), collapse = " + "))))

  descadd <- forms1[!chosen]
  dropadd <- rep("add", sum(!chosen))

  #build a formula list dropping elements
  if (any(chosen)) {
    if (sum(chosen)==1) {
      formsdrop <- list(n ~ 1)
    } else {
      formsdrop <- 
        lapply(seq(sum(chosen)), 
           function(i) 
             as.formula(paste0("n ~ ", paste(forms1[chosen][-i], collapse = " + ") )))
    }
    descdrop <- forms1[chosen]
    forms <- c(formsadd, formsdrop)
    desc <- c(descadd, descdrop)
    dropadd <- c(dropadd, rep("drop", sum(chosen)))
  } else {
    forms <- formsadd
    desc <- descadd
  }

  mods <- lapply(forms, efp, data = data, pass = data$pass)

  if (all(!chosen)) {
    m0 <- efp(X ~ 1, data = data, pass = data$pass)
  } else {
    m0 <- efp(as.formula(paste0("n ~ ", paste(forms1[chosen], collapse = " + "))), data = data, pass = data$pass)    
  }
  
  tab <- summaryMods(mods, m0 = m0, order = FALSE, fn = fn, phi = phi, nobs = nobs)[,-1]
  tab <- cbind(what = desc, step = dropadd,  tab)
  tab[order(tab $ Dic),]
}
  

run1Selection <- function(forms1, data = ef, fn = phiBIC, 
                          phi = 2.981, nobs = 2838, start = NULL,
                          tol = 1e-3) {
  chosen <- rep(FALSE, length(forms1))
  if (!is.null(start)) chosen <- start

  out <- runModels(chosen, forms1, data = data, fn = fn, phi = phi, nobs = nobs)
  print(head(out, 10), digits = 3)
  if ( out $ Dic[1] < -tol) {
    # then model is an improvement
    which <- which(forms1 %in% out $ what[1])
    chosen[which] <- !chosen[which]
  }
  list(tab = out, chosen = chosen)
}


runSelection <- function(forms1, data = ef, fn = phiBIC, phi = 2.981, nobs = 2838, start = NULL) {
  chosen <- rep(FALSE, length(forms1))
  if (!is.null(start)) chosen <- start
  out <- list(Dic = -1)
  tol <- 0
  outlist <- list() # grow this - bad!
  i <- 0
  while(out $ Dic[1] < tol & !all(chosen)) {
    out <- runModels(chosen, forms1, data = data, fn = fn, phi = phi, nobs = nobs)
    print(head(out, 10), digits = 3)
    outlist[[i <- i + 1]] <- out
    if ( out $ Dic[1] < 0) {
      # then model is an improvement
      which <- which(forms1 %in% out $ what[1])
      chosen[which] <- !chosen[which]
      cat("\n ----------------------------------------------- \n\t", paste(forms1[chosen], collapse = " + "), "\n ----------------------------------------------- \n\n")
    } 
  }

  if (all(!chosen)) {
    mod <- efp(n ~ 1, data = data, pass = data$pass)
  } else {
    mod <- efp(as.formula(paste0("n ~ ", paste(forms1[chosen], collapse = " + "))), data = data, pass = data$pass)
  }

  list(history = outlist, chosen = chosen, mod = mod)
}


