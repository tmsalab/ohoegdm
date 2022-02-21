# Data Analysis for PISA 12 US Vignette ---- 

## Load data ----
y = edmdata::items_ordered_pisa12_us_vignette

## Configure Q and BETA parameters ----
K = 3        # Number of Attributes
order = K    # Maximum number of interactions included in Theta/Beta. Default is K.
M = 4        # Number of Response Categories

## Chain Properties ----
burnin = 20000 
chain_length = 80000

## Spike-slab callibration ---- 
l0s = c(1, rep(100, 2 ^ K - 1))
l1s = c(1, rep(0.5, 2 ^ K - 1))

## Fit model ---- 
fit_model =
  ohoegdm::ohoegdm(
    y = y,
    k = K,
    m = M,
    order = order,
    sd_mh = .04,
    burnin = burnin,
    chain_length = chain_length,
    l0 = l0s,
    l1 = l1s
  )


## Analyze model results ---- 

### Check MH acceptance
fit_model$recoovery$MHsum

### Disabled on release, please uncomment and recompile in src/ if needed
### PPP results 
# item_mean_ppp = fit_model$chain$itemmeanppp
# item_mean_cov_ppp = fit_model$chain$itemcovppp
# summary(c(fit_model$item_mean_ppp))

# upelem = item_mean_cov_ppp[upper.tri(item_mean_cov_ppp, diag = T)]
# 1 - sum(1 * (upelem < .025 | upelem > .975)) / length(upelem)

### Absolute fit measures
# mean(fit_model$chain$srmr)

# (seek help,exhaust individual knowledge,act quickly)
# Threshold the estimate Q matrix into a dichotomous 0/1 entry. 
Q_est = 1*(fit_model$estimates$estimates$QS > .5)

# View parameter estimates
mbeta = round(fit_model$estimates$betas$mean, 2)
mkappas= round(fit_model$estimates$kappas$mean, 3)
mtaus = round(fit_model$estimates$taus$mean, 2)
mlambdas = round(apply(fit_model$chain$lambdas, 1, mean), 2) # Consider rotating.
mthetas = fit_model$estimates$mtheta
sdmtheta = sd(mthetas)
