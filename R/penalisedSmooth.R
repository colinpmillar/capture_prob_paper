
library(CLmodel)


#' Estimate capture probabilites from electrofishing data
#'
#' This function uses the marginal likelihood of capture probabilities
#' to estimate model parameters 
#' 
#'
#' @param formula a formula object
#' @param data a data.frame containing all relavent info
#' @return glm type object
#' @export
#' @examples
#' # none yet
efp <- function(formula, data = NULL, passes = NULL, verbose=TRUE, init = "0", hessian = FALSE) {

  if (!exists("stanmod")) {
    message("Building optimiser for first use...")
    stanmod <- rstan::stan_model(model_code = "
      data {
        int<lower=0> N; // number of observations
        int<lower=0> K; // number of parameters
        real S[N]; // the number of fishing passes
        real R[N]; // Zippins R (see seber p 312, eq 7.22)
        real T[N]; // total catches
        matrix[N,K] A;
      }
      parameters {
        vector[K] alpha;
      } 
      model {
        vector[N] expeta;
        expeta <- exp(A * alpha);
        for (i in 1:N) {
          real p;
          p <- expeta[i]/(1.0 + expeta[i]);
          increment_log_prob(T[i] * log(p));
          increment_log_prob(T[i] * R[i] * log(1-p));
          increment_log_prob(-T[i] * log(1 - (1-p)^S[i]) );
        }
      }")
    assign("stanmod", stanmod, .GlobalEnv)
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


  standat <- 
    list(N = nrow(Gfit), K = ncol(Gfit), 
         S = data0 $ S, T = data0 $ T, R = with(data0, S - 1 - Z),
         A = Gfit)
  if (!verbose) {
    tmp <- 
      capture.output(
        opt <- rstan::optimizing(stanmod, data = standat, algorith = "BFGS", hessian = hessian, verbose = verbose, init = init)
      )
  } else {
    opt <- rstan::optimizing(stanmod, data = standat, algorith = "BFGS", hessian = hessian, verbose = verbose, init = init)
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

