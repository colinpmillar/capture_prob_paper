
test <- function(formula = Z ~ s(year) + s(doy), lambda, 
                 data = wk, passes = "Runs", 
                 hessian = TRUE, verbose = TRUE, init = "0") {

  
  if (!exists("stanmod2")) {
    message("Building optimiser for first use...")
    stanmod2 <- rstan::stan_model(model_code = "
      data {
        int<lower=0> N; // number of observations
        int<lower=0> K; // number of parameters
        real S[N]; // the number of fishing passes
        real R[N]; // Zippins R (see seber p 312, eq 7.22)
        real T[N]; // total catches
        matrix[N,K] A; // model matrix
        matrix[K,K] Q; // penalty matrix
      }
      parameters {
        vector[K] alpha;
      } 
      model {
        vector[N] expeta;
        vector[N] p;
        expeta <- exp(A * alpha);
        p <- expeta ./ (1.0 + expeta);
        for (i in 1:N) {
          increment_log_prob(T[i] * log(p[i]));
          increment_log_prob(T[i] * R[i] * log(1-p[i]));
          increment_log_prob(-T[i] * log(1 - (1-p[i])^S[i]) );
        }
        increment_log_prob(-1.0 * quad_form(Q, alpha));
      }")
    assign("stanmod2", stanmod2, .GlobalEnv)
  }

  if (is.null(data)) stop("must supply data")
  data0 <- subset(data, T > 0)

  if (is.null(passes)) stop("must supply the number of fishing runs")
  data0 $ S <- data0[[passes]]

  # set up model
  if (nrow(data0) == 1) {
    G <- matrix(1, 1, 1)
  } else {
    Gsetup <- gam(formula, data = data0, fit = FALSE)
    G <- Gsetup $ X
  }

  # remove redundant / collinear parameters
  qr.G <- qr(G)
  rank.deficient <- qr.G $ pivot[abs(diag(qr.G $ qr)) < 1e-7]
  if (length(rank.deficient)) {
    droppar <- paste(colnames(G)[rank.deficient], collapse = "\n\t")
    warning("*** Model has ", length(rank.deficient)," too many parameter(s)!!\n    i will remove the redundant ones:\n\t", droppar, call. = FALSE)
    Gfit <- G[,-rank.deficient]
  } else {
    Gfit <- G
  }

  # build up penalty matrix
  ## for now make it the zero matrix for unpenalised
  Q <- matrix(0, ncol(Gfit), ncol(Gfit))
  nsmooth <- length(Gsetup $ S)
  if (length(lambda) != nsmooth) stop("lambda needs to be", nsmooth, "in length!")
  for (i in 1:nsmooth) {
    ind <- Gsetup $ smooth[[i]] $ first.para:Gsetup $ smooth[[i]] $ last.para
    Q[ind,ind] <- lambda[i] * Gsetup $ S[[i]]
  }  
  
  standat <- 
    list(N = nrow(Gfit), K = ncol(Gfit), 
         S = data0 $ S, T = data0 $ T, R = with(data0, S - 1 - Z),
         A = Gfit, 
         Q = Q)

  if (!verbose) {
    tmp <- 
      capture.output(
        opt <- rstan::optimizing(stanmod2, data = standat, algorith = "BFGS", hessian = hessian, verbose = verbose, init = init)
      )
  } else {
    opt <- rstan::optimizing(stanmod2, data = standat, algorith = "BFGS", hessian = hessian, verbose = verbose, init = init)
  } 

  opt $ formula <- formula # for printing and summary
  opt $ llik <- opt $ value
  opt $ terms <- Gsetup $ terms
  opt $ call <- match.call()
  opt $ aic <- -2 * opt $ llik + 2 * ncol(Gfit)
  opt $ G <- G
  opt $ Gfit <- Gfit
  opt $ coefficients <- opt $ par
  names(opt $ coefficients) <- colnames(Gfit)
  opt $ df.null <- nrow(G)
  opt $ df.residual <- nrow(G) - ncol(Gfit)
  opt $ rank <- ncol(Gfit)
  opt $ fitted <- p <- transpar(opt $ par, Gfit)
  opt $ residuals <- rep(0, nrow(data0))
  opt $ null.deviance <- NA
  opt $ deviance <- NA 
  opt $ family <- binomial()
  opt $ Vb <- if (hessian) try(solve(-1 * opt $ hessian)) else NULL
  opt $ Gsetup <- Gsetup

  # calculate effetive degrees of freedom
  

  # get a gam container
  # g1 <- gam(G = Gsetup)
  # g1 $ coefficients[] <- opt $ par       
  # g1 $ Vp[] <- opt $ Vb
  # g1 $ family <- binomial()
  # X <- predict(g1, type = "lpmatrix")
  # g1 $ linear.predictors <-  c(X %*% g1 $ coef)
  # g1 $ fitted.values <- c(1/(1 + exp(-g1 $ linear.predictors)))
  # g1 $ aic <- opt $ aic

  class(opt) <- c("efp", "glm", "lm")
  opt
}


# predict p for data
predictX <- function(mod, newdata) {
  g1 <- gam(G = mod $ Gsetup)
  qr.G <- qr(mod $ G)
  rank.deficient <- qr.G $ pivot[abs(diag(qr.G $ qr)) < 1e-7]
  whichkeep <- -rank.deficient
  if (!length(whichkeep)) whichkeep <- 1:length(mod $ coefficients) 
  names(g1 $ coefficients[-1 * whichkeep])
  g1 $ coefficients[] <- 0
  g1 $ coefficients[whichkeep] <- mod $ coefficients       
  g1 $ Vp[] <- 0
  diag(g1 $ Vp[]) <- 1e-5
  g1 $ Vp[whichkeep, whichkeep] <- mod $ Vb
  g1 $ family <- binomial()

  predict(g1, type = "lpmatrix", newdata = newdata)
}





