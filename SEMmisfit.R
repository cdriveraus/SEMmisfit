# ================================================================
# Simulation Study for SEM Global Model Fit Tests with HFI
# ================================================================
#
# This script performs a simulation study to evaluate the performance
# of traditional SEM chi-square tests versus nonlinear alternatives.
# The study includes varying sample sizes and multiple misfit types:

required_packages <- c(
  "MASS", # For mvrnorm function
  "OpenMx", # For SEM modeling
  "data.table", # For data manipulation
  "foreach", # For parallel processing
  "doParallel", # For parallel processing
  "ggplot2", # For data visualization
  "moments", # For skewness and kurtosis
  "goftest", # For goodness-of-fit tests
  'energy', # For mvnorm.etest
  'infotheo' # For mutinformation
)

installed_packages <- rownames(installed.packages())

for(pkg in required_packages){ # Install required packages if not already installed
  if(!pkg %in% installed_packages){
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}



# Detect the number of available cores
num_cores <- detectCores() - 1  # Reserve one core for the OS

# Register parallel backend
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# ----------------------------
#  Define Simulation Parameters
# ----------------------------

niter <- 1000 # Number of iterations per condition

simconditions <- data.table(expand.grid(
  n = c(100, 500),
  nvars=c(3,5),
  misfit_type = c('none', 'classical', 'groupvar','groupmean', 'interaction','nonnormal'),
  stringsAsFactors = FALSE
))

simresults <- data.table()  # Initialize an empty data.table to store results


# ----------------------------
#  Define the Simulation Function
# ----------------------------

# Function to perform a single simulation iteration
simulate_iteration <- function(cond, iter){
  
  # Extract condition parameters
  #cond <- simconditions[1,]
  n <- cond$n
  misfit_type <- cond$misfit_type
  nvars <- cond$nvars
  
  # -----------------------------------
  # Data Generation
  # -----------------------------------

  data <- data.table(mvrnorm(n, mu = rep(0, nvars), Sigma = diag(1,nvars)))
  
  
  data[,V3:=  0 * V1 + 
      0.3 * V2 + 
      ifelse('interaction' %in% misfit_type, 1, 0) * (V1 * V2) + # Introduce nonlinearity (if needed) by an interaction term
      rnorm(n, mean = 0, sd = 1)
  ]
  
  if('groupvar' %in% misfit_type) data[1:floor(n/2),c('V2','V3') ] <- data[1:(n/2),c('V2','V3'),with=F ] * 2
  if('groupmean' %in% misfit_type) data[1:floor(n/2),c('V2','V3') ] <- data[1:(n/2),c('V2','V3'),with=F ] + 2
  
  if('nonnormal' %in% misfit_type) data[,V1:=log1p(exp(V1))]
  
  # -----------------------------------
  #  Define the SEM Model
  # -----------------------------------
  
  # Define manifest variables
  manifests <- paste0('V',1:nvars)
  
  # Specify the SEM path model with fixed covariance between X1 and X2
  model <- mxModel("PathModel", #define free covariances except between x1 and x2
    type = "RAM",
    manifestVars = manifests,
    latentVars = NULL,
    mxPath(from = "V1", to = "V3", arrows = 2, free = TRUE, values = 0.5, labels = "b1"),
    mxPath(from = "V2", to = "V3", arrows = 2, #if there is classical covariance misfit, fix covariance between X2 and X3 to incorrect value
      free = ifelse('classical' %in% misfit_type,FALSE,TRUE), values = 0.8, labels = "b2"),
    mxPath(from = manifests, arrows = 2, free = T, values = 1),
    mxPath(from = 'one', to = manifests),
    mxFitFunctionML(rowDiagnostics = TRUE),
    mxData(observed = data, type = "raw")
  )
  
  # -----------------------------------
  #  Fit the SEM Model
  # -----------------------------------
  
  # Fit the model with error handling
  fit <- tryCatch({
    suppressMessages(mxRun(model))
  }, error = function(e) {
    # Return NA for p-values if model fails to run
    return('fit error')
  })
  
  # Check if model fitting was successful
  if(inherits(fit, "MxModel") == FALSE){
    return('fit error check')
  }
  
  # Check for model convergence
  if(fit$output$status$code != 0){
    return('non convergence')
  }
  
  # -----------------------------------
  #  Compare with Reference Model
  # -----------------------------------
  
  # Generate reference (saturated) models
  refs <- tryCatch({
    suppressMessages(mxRefModels(fit, run = TRUE))
  }, error = function(e){
    return('ref model error')
  })
  
  
  # Extract and store the SEM chi-square test p-value
  SEM_p <-  mxCompare(refs[[1]], fit)$p
  
  # -----------------------------------
  #  Extract Residuals and Perform Misfit Tests
  # -----------------------------------
  
  # Extract residuals
  expMean = c(attributes(fit$output$algebras$PathModel.fitfunction)$expMean)
  expCov <- attributes(fit$output$algebras$PathModel.fitfunction)$expCov
  residuals <- t(apply(data, 1, function(row) row - expMean))
  
  # Extract row-wise sum of squared standardized residuals
  stdresidualsrowsum <- attributes(fit$output$algebras$PathModel.fitfunction)$rowDist
  
  # Function to compute conditional predictions and standardized residuals
  conditional_residuals <- {
    conditional_residuals <- matrix(NA, nrow = nrow(data), ncol = nvars)
    for (j in 1:nvars) {
      # Partition covariance matrix
      Sigma_jj <- expCov[j, j, drop=TRUE]  # Extract variance of X_j (scalar)
      Sigma_j_rest <- expCov[j, -j, drop=FALSE]  # Extract row vector (1 × (p-1))
      Sigma_rest_rest <- expCov[-j, -j, drop=FALSE]  # Extract covariance matrix ((p-1) × (p-1))
      
      # Invert the covariance matrix of remaining variables
      Sigma_inv <- solve(Sigma_rest_rest)  # Inverse of (p-1) × (p-1) matrix
      
      # Compute conditional variance (scalar)
      cond_var <- Sigma_jj - Sigma_j_rest %*% Sigma_inv %*% t(Sigma_j_rest)  
      
      # Compute conditional mean of X_j given all other variables
      X_rest <- data[, colnames(data)[-j], with=FALSE]  # Extract all variables except X_j
      mu_j <- expMean[j] + as.vector(Sigma_j_rest %*% Sigma_inv %*% t(X_rest - expMean[-j]))  # Convert to vector
      
      # Compute standardized residuals
      conditional_residuals[, j] <- c((data[[ colnames(data)[j]]] - mu_j) / c(sqrt(cond_var)))
    }
    conditional_residuals
  }

  
  ##manually compute squared standardized residuals rowsum using expectations
  # stdresiduals <- apply(residuals,1,function(x) x %*% solve(expCov) %*% x)
  
  df <- length(manifests)  # Degrees of freedom for ks chi-square test (number of observed variables)
  
  # Perform tests of std residuals against chi-square distribution
  AD_p <- ad.test(stdresidualsrowsum, 'pchisq', df = df)$p.value
  
  CVM_p <-cvm.test(stdresidualsrowsum, 'pchisq', df = df)$p.value
  
  KS_p <- ks.test(stdresidualsrowsum, 'pchisq', df = df)$p.value
  
  univarAD_p <- {
    adtest <- c()
    for(vari in 1:nvars){
      adtest <- c(adtest, ad.test(residuals[,vari], 'pnorm',expMean[vari],sqrt(expCov[vari,vari]))$p.value)
    }
    min(adtest)*nvars #bonferroni correction
  }
  
  
  MVnorm_p <- mvnorm.etest(residuals,1000)$p.value
  

  
  # # Levene’s test for heteroscedasticity in residual variance
  # Levene_p <- {
  #   levenetest <- c()
  #   for(vari in 1:nvars){
  #     # Grouping continuous predictors into quantiles
  #     group_var <- cut(residuals[, vari], breaks = 3)  
  #     levenetest <- c(levenetest, leveneTest(stdresidualsrowsum ~ group_var)$"Pr(>F)"[1])
  #   }
  #   min(levenetest) * nvars # Bonferroni correction
  # }
  
  # Mutual Information test for nonlinear dependencies
  MutInfo_p <- {
    mi_results <- c()
    p_values <- c()
    for (vari in 1:nvars) {
      # Compute observed MI
      observed_mi <- mutinformation(discretize(apply(conditional_residuals[, -vari,drop=F],1,function(x) sum(x^2))), discretize(conditional_residuals[, vari]))
      # Generate null distribution by permuting R
      permuted_mis <- replicate(1000, {
        mutinformation(discretize(sample(stdresidualsrowsum)), discretize(conditional_residuals[, vari]))
      })
      # Compute p-value: Proportion of permuted MIs greater than or equal to observed MI
      p_values <- c(p_values, mean(permuted_mis >= observed_mi))
    }
    # Apply Bonferroni correction for multiple comparisons
    corrected_p_value <- min(p_values) * nvars
  }
  
  
  # -----------------------------------
  #  Compute Heteroscedasticity Fit Index (HFI)
  # -----------------------------------
  
  compute_HFI <- function(residuals) {
    n <- length(residuals)
    
    # Compute sums of squares and fourth powers
    sum_e2 <- sum(residuals^2)   # Sum of squared residuals (scalar)
    sum_e4 <- sum(residuals^4)   # Sum of residuals to the fourth power (scalar)
    
    # Avoid division by zero (small or zero residual variance)
    if (sum_e2 == 0) {
      return(1)  # Perfect fit: no residual variance
    }
    
    # Compute heteroscedasticity measure (h_het)
    h_het <- sqrt(n / 24) * (
      ((1/n) * sum_e4) / 
        ((1/(n) * sum_e2)^2) - 
        3)
    
    # Adjust h_het to compute HFI
    if (h_het <= 0) {
      return(1)  # Homoscedastic case
    } else {
      a <- 0.032  # Scaling factor
      return(1 / (a * h_het + 1))
    }
  }
  
  HFI=compute_HFI(residuals)
  
  # Return p-values and HFI
  return(data.table(SEM_p = SEM_p, KS_p = KS_p, HFI = HFI, AD_p = AD_p,CVM_p=CVM_p,univarAD_p=univarAD_p, MVnorm_p=MVnorm_p,MutInfo_p=MutInfo_p))
}

# ----------------------------
#  Run Simulations in Parallel
# ----------------------------

# Use foreach to iterate over each simulation condition
# Outer loop iterates over conditions
# Inner loop iterates over iterations within each condition
simresultslist <- foreach(cond_idx = 1:nrow(simconditions), #.combine = rbind, 
  .packages = c("OpenMx", "MASS", "data.table", "moments", "goftest",'energy', 'car', 'infotheo')) %:%
  foreach(iter = 1:niter, .packages = c("OpenMx", "MASS", "data.table", "moments")) %dopar% {
    
    # Extract current condition
    cond <- simconditions[cond_idx]
    
    # Perform simulation iteration
    data.table(condition = cond_idx,
      n = cond$n,
      nvars=cond$nvars,
      misfit_type = cond$misfit_type,
      simulate_iteration(cond, iter)
    )
  }

simresultslist <- unlist(simresultslist, recursive = FALSE)
outlength = unlist(lapply(simresultslist, length))
if(any(outlength != max(outlength))){
  warning('errors occurred for some iterations')
  print(simresultslist[which(outlength != max(outlength))])
  simresultslist <- simresultslist[which(outlength == max(outlength))]
}
simresults <- rbindlist(simresultslist)

# Shut down the parallel cluster
stopCluster(cl)
registerDoSEQ()

# ----------------------------
# Analyze and Summarize Results
# ----------------------------

# Summary statistics for each condition
summary_results <- simresults[, .(
  SEM = round(mean(SEM_p <= 0.05, na.rm = TRUE), 3),
  KS = round(mean(KS_p <= 0.05, na.rm = TRUE), 3),
  univarAD = round(mean(univarAD_p <= 0.05, na.rm = TRUE), 3),
  AD = round(mean(AD_p <= 0.05, na.rm = TRUE), 3),
  CVM = round(mean(CVM_p <= 0.05, na.rm = TRUE), 3),
  MVnorm = round(mean(MVnorm_p <= 0.05, na.rm = TRUE), 3),
  MutInfo = round(mean(MutInfo_p <= 0.05, na.rm = TRUE), 3),
  HFI = round(mean(HFI <= 0.95, na.rm = TRUE), 3)  
), by = .(condition, n, nvars, misfit_type)]

# Print summary results
print(summary_results)

## Optionally, save the results to CSV files
# fwrite(simresults, "SEM_Simulation_Results.csv")
# fwrite(summary_results, "SEM_Simulation_Summary.csv")

# ----------------------------
#  Visualization 
# ----------------------------

summary_long <- melt(
  summary_results,
  id.vars = c("condition", colnames(simconditions)),
  measure.vars = colnames(summary_results)[!colnames(summary_results) %in% c("condition",'n','misfit_type','nvars')],
  variable.name = "misfit_metric",
  value.name = "misfit_rate"
)


summary_long[,`Misfit Type`:=misfit_type]
# Plot misfit rates for SEM chi-square test, Raw KS test, and HFI across conditions
p <- ggplot(
  summary_long[!misfit_metric %in% c('KS','CVM','Levene'),],
  aes(x = factor(n), y = misfit_rate, fill = misfit_metric)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  facet_grid(rows = vars(`Misfit Type`),cols=vars(nvars),labeller = label_both) +
  labs(
    x = "Sample Size (n)",
    y = "Proportion of p-values ≤ 0.05",
    fill = "Misfit Metric"
  ) +
  theme_bw() +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
  theme(legend.position = "bottom") 

print(p)
pdf('sim.pdf', width = 10, height = 12)
print(p)
dev.off()

# ================================================================
# End of Simulation Study Script
# ================================================================
