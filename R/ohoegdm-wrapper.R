# Construct the OHOEGDM object
new_ohoegdm_model = function(model_mcmc,
                             details,
                             estimates = NULL,
                             recovery = NULL) {

    estimate_chain = lapply(model_mcmc, summarize_model)
    estimates = c(estimate_chain, estimates)
    chain = model_mcmc[!(names(model_mcmc) %in% estimates)]
    structure(
        list(
            estimates = estimates, # Iterates over each parameter and obtains summary information
            chain = chain,
            details = details,
            recovery = recovery
        ),
        class = c("ohoegdm", "egdm")
    )
}


#' Ordinal Higher-Order General Diagnostic Model under the
#' Exploratory Framework (OHOEGDM)
#'
#' Performs the Gibbs sampling routine for an ordinal higher-order EGDM.
#'
#' @param y                 Ordinal Item Matrix
#' @param k                 Dimension to estimate for Q matrix
#' @param m                 Number of Item Categories. Default is `2` matching the binary case.
#' @param order             Highest interaction order to consider. Default model-specified `k`.
#' @param sd_mh             Metropolis-Hastings standard deviation tuning parameter. 
#' @param burnin            Amount of Draws to Burn
#' @param chain_length      Number of Iterations for chain.
#' @param l0                Spike parameter. Default 1 for intercept and 100
#'                          coefficients
#' @param l1                Slab parameter. Default 1 for all values. 
#' @param m0,bq             Additional tuning parameters.
#' @return
#'
#' A `ohoegdm` object containing four named lists:
#'
#' - **`estimates`**: Averaged chain iterations
#'    - `thetas`: Average theta coefficients
#'    - `betas`: Average beta coefficients
#'    - `deltas`: Average activeness of coefficients
#'    - `classes`: Average class membership
#'    - `m2lls`: Average negative two times log-likelihood
#'    - `omegas`: Average omega
#'    - `kappas` : Average category threshold parameter
#'    - `taus`: Average \eqn{K}-vectors of factor intercept
#'    - `lambdas`: Average \eqn{K}-vectors of factor loadings
#'    - `guessing`: Average guessing item parameter
#'    - `slipping`: Average slipping item parameter
#'    - `QS`: Average activeness of Q matrix entries
#' - **`chain`**: Chain iterations from the underlying _C++_ rountine.
#'    - `thetas`: Theta coefficients iterations
#'    - `betas`:  Beta coefficients iterations
#'    - `deltas`: Activeness of coefficients iterations
#'    - `classes`:  Class membership iterations
#'    - `m2lls`: Negative two times log-likelihood iterations
#'    - `omegas`:  Omega iterations
#'    - `kappas` : Category threshold parameter iterations
#'    - `taus`: \eqn{K}-vector of factor intercept iterations
#'    - `lambdas`: \eqn{K}-vector of factor loadings iterations
#'    - `guessing`: Guessing item parameter iterations
#'    - `slipping`: Slipping item parameter iterations
#' - **`details`**: Properties used to estimate the model
#'    - `n`: Number of Subjects
#'    - `j`: Number of Items
#'    - `k`: Number of Traits
#'    - `m`: Number of Item Categories.
#'    - `order`: Highest interaction order to consider. Default model-specified `k`.
#'    - `sd_mh`: Metropolis-Hastings standard deviation tuning parameter. 
#'    - `l0`: Spike parameter
#'    - `l1`: Slab parameter
#'    - `m0`, `bq`: Additional tuning parameters
#'    - `burnin`: Number of Iterations to discard
#'    - `chain_length`: Number of Iterations to keep
#'    - `runtime`: Elapsed time algorithm run time in the _C++_ code.
#' - **`recovery`**: Assess recovery metrics under a simulation study.
#'    - `Q_item_encoded`: Per-iteration item encodings from Q matrix.
#'    - `MHsum`: Average acceptance from metropolis hastings sampler
#' 
#' @details
#'
#' The **`estimates`** list contains the mean information from the sampling
#' procedure. Meanwhile, the **`chain`** list contains full MCMC values. Moreover,
#' the **`details`** list provides information regarding the estimation call.
#' Lastly, the **`recovery`** list stores values that can be used when
#' assessing the method under a simulation study.
#'
#' @rdname ohoegdm
#' @export
#'
#' @examples
#' # Simulation Study
#' if (requireNamespace("edmdata", quietly = TRUE)) {
#' # Q and Beta Design ----
#' 
#' # Obtain the full K3 Q matrix from edmdata
#' data("qmatrix_oracle_k3_j20", package = "edmdata")
#' Q_full = qmatrix_oracle_k3_j20
#' 
#' # Retain only a subset of the original Q matrix
#' removal_idx = -c(3, 5, 9, 12, 15, 18, 19, 20)
#' Q = Q_full[removal_idx, ]
#' 
#' # Construct the beta matrix by-hand
#' beta = matrix(0, 20, ncol = 8)
#' 
#' # Intercept
#' beta[, 1] = 1
#' 
#' # Main effects
#' beta[1:3, 2] = 1.5
#' beta[4:6, 3] = 1.5
#' beta[7:9, 5] = 1.5
#' 
#' # Setup two-way effects
#' beta[10, c(2, 3)] = 1
#' beta[11, c(3, 4)] = 1
#' 
#' beta[12, c(2, 5)] = 1
#' beta[13, c(2, 5)] = 1
#' beta[14, c(2, 6)] = 1
#' 
#' beta[15, c(3, 5)] = 1
#' beta[16, c(3, 5)] = 1
#' beta[17, c(3, 7)] = 1
#' 
#' # Setup three-way effects
#' beta[18:20, c(2, 3, 5)] = 0.75
#' 
#' # Decrease the number of Beta rows
#' beta = beta[removal_idx,]
#' 
#' # Construct additional parameters for data simulation
#' Kappa = matrix(c(0, 1, 2), nrow = 20, ncol = 3, byrow =TRUE) # mkappa
#' lambda = c(0.25, 1.5, -1.25) # mlambdas
#' tau = c(0, -0.5, 0.5) # mtaus
#' 
#' 
#' # Simulation conditions ---- 
#' N = 100        # Number of Observations
#' J = nrow(beta) # Number of Items
#' M = 4          # Number of Response Categories
#' Malpha = 2     # Number of Classes
#' K = ncol(Q)    # Number of Attributes
#' order = K      # Highest interaction to consider
#' sdmtheta = 1   # Standard deviation for theta values
#' 
#' # Simulate data ---- 
#' 
#' # Generate theta values
#' theta = rnorm(N, sd = sdmtheta)
#' 
#' # Generate alphas 
#' Zs = matrix(1, N, 1) %*% tau + 
#'      matrix(theta, N, 1) %*% lambda + 
#'      matrix(rnorm(N * K), N, K)
#' Alphas = 1 * (Zs > 0)
#' 
#' vv = gen_bijectionvector(K, Malpha)
#' CLs = Alphas %*% vv
#' Atab = GenerateAtable(Malpha ^ K, K, Malpha, order)$Atable
#' 
#' # Simulate item-level data
#' Ysim = sim_slcm(N, J, M, Malpha ^ K, CLs, Atab, beta, Kappa)
#' 
#' # Establish chain properties 
#' # Standard Deviation of MH. Set depending on sample size.
#' # If sample size is:
#' #  - small, allow for larger standard deviation
#' #  - large, allow for smaller standard deviation.
#' sd_mh = .4 
#' burnin = 50        # Set for demonstration purposes, increase to at least 5,000 in practice.
#' chain_length = 100 # Set for demonstration purposes, increase to at least 40,000 in practice.
#' 
#' # Setup spike-slab parameters
#' l0s = c(1, rep(100, Malpha ^ K - 1))
#' l1s = c(1, rep(1, Malpha ^ K - 1))
#' 
#' my_model = ohoegdm::ohoegdm(
#'   y = Ysim,
#'   k = K,
#'   m = M,
#'   order = order,
#'   l0 = l0s,
#'   l1 = l1s,
#'   m0 = 0,
#'   bq = 1,
#'   sd_mh = sd_mh,
#'   burnin = burnin,
#'   chain_length = chain_length
#' )
#' }
ohoegdm = function(y,
                 k,
                 m = 2,
                 order = k,
                 sd_mh = 0.4,
                 burnin = 1000L,
                 chain_length = 10000L,
                 l0 = c(1, rep(100, sum(choose(k, seq_len(order))))),
                 l1 = c(1, rep(1, sum(choose(k, seq_len(order))))),
                 m0 = 0,
                 bq = 1) {

    # Perform some quality checks
    stopifnot(is.matrix(y))
    stopifnot(chain_length > burnin)

    # parameter_length = sum(choose(k, seq_len(order))) + 1
    # stopifnot(length(l0) == parameter_length)
    # stopifnot(length(psi_invj) == parameter_length)
    stopifnot(order <= k & order >= 1)

    # Launch routine and time it
    timing = system.time({
        model_mcmc <- ohoegdm_cpp(
            Y = y, K = k, M = m, order = order,
            l0 = l0, l1 = l1,
            m0 = m0, bq = bq,
            sdMH = sd_mh,
            burnin = burnin, chain_length = chain_length)
    })

    
    # Package object
    new_ohoegdm_model(
        model_mcmc$chain,
        estimates = model_mcmc$estimates,
        details = list(
            n = nrow(y),
            j = ncol(y),
            k = k,
            m = m,
            order = order,
            sd_mh = sd_mh,
            l0 = l0,
            l1 = l1,
            m0 = m0,
            bq = bq,
            burnin = burnin,
            chain_length = chain_length,
            runtime = timing[["elapsed"]]
        ),
        recovery = model_mcmc$recovery
    )
}


