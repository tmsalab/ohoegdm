# Construct the EXMS object
new_ohoegdm_model = function(model_mcmc,
                             details,
                             estimates = NULL) {

    estimate_chain = lapply(model_mcmc, summarize_model)
    estimates = c(estimate_chain, estimates)
    chain = model_mcmc[!(names(model_mcmc) %in% estimates)]
    structure(
        list(
            estimates = estimates, # Iterates over each parameter and obtains summary information
            chain = chain,
            details = details
        ),
        class = c("ohoegdm", "egdm")
    )
}


#' Ordinal Higher-Order General Diagnostic Model under the
#' Exploratory Framework (OHOEGDM)
#'
#' Performs the Gibbs sampling routine for a higher-order EGDM.
#'
#' @param y                 Item Matrix
#' @param k                 Dimension to estimate for Q matrix
#' @param m                 Number of Item Categories. Default is `2` matching the binary case.
#' @param order             Highest interaction order to consider. Default model-specified `k`.
#' @param sdMH              Metropolis-Hastings standard deviation tuning parameter. 
#' @param burnin            Amount of Draws to Burn
#' @param chain_length      Number of Iterations for chain.
#' @param l0,l1,m0,bq Additional tuning parameters.
#' @return
#'
#' A `egdm_delta` object containing three named lists:
#'
#' **`estimates`**
#'
#' - `thetas`: Average theta coefficients
#' - `betas`: Average beta coefficients
#' - `class`: Average class membership
#' - `pi`: Average attribute class probability.
#' - `omega`: Average omega
#'
#' **`chain`**
#'
#' - `thetas`: theta coefficients iterations
#' - `betas`:  beta coefficients iterations
#' - `class`:  class membership iterations
#' - `omega`:  omega iterations
#'
#' **`details`**
#'
#' - `n`: Number of Subjects
#' - `j`: Number of Items
#' - `k`: Number of Traits
#' - `order`: Highest interaction order to consider. Default model-specified `k`.
#' - `sdMH`: Metropolis-Hastings standard deviation tuning parameter. 
#' - `l0`: Spike parameter
#' - `l1`: Slab parameter
#' - `m0`, `bq`: Additional tuning parameters
#' - `burnin`: Number of Iterations to discard
#' - `chain_length`: Number of Iterations to keep
#'
#' @details
#'
#' The **`estimates`** list contains the mean information from the sampling
#' procedure. Meanwhile, the **`chain`** list contains full MCMC values. Lastly,
#' the **`details`** list provides information regarding the estimation call.
#'
#' @rdname ohoegdm
#' @export
#'
#' @examples
#' \dontrun{
#' ## Configuration
#' N = 100  # Number of Subjects
#' J = 20   # Number of Items
#' K = 3    # Number of Attributes
#'
#' ## TBA
#' }
ohoegdm = function(y,
                 k,
                 m = 2,
                 order = k,
                 sdMH = 0.4,
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
            sdMH = sdMH,
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
            sdMH = sdMH,
            l0 = l0,
            l1 = l1,
            m0 = m0,
            bq = bq,
            burnin = burnin,
            chain_length = chain_length,
            runtime = timing[["elapsed"]]
        )
    )
}


