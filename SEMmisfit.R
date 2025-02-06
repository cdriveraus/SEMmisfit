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
  'infotheo', # For mutinformation
  'dHSIC' # For distance‐based Hilbert–Schmidt independence criterion
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
nboot <- 1000 # Number of bootstrap samples for permutation tests

simconditions <- data.table(expand.grid(
  n = c(100,500),
  nvars=c(4,8),
  misfit_type = c('none', 'classical', 'interaction','nonnormal', 'groupvar','groupmean'),
  stringsAsFactors = FALSE
))

simresults <- data.table()  # Initialize an empty data.table to store results

# Function to generate signed chi-square distributed variables / cumulative distribution
generate_signed_chisq <- function(n, df = 3) {
  X <- matrix(rnorm(n * df), ncol = df)  # Generate standard normal variables
  rowSums(sign(X) * X^2)  # Apply transformation
}

#function to combine multiple p values into one
combine_pvalues <- function(p_values, method = c("Fisher")) {
  if(any(is.na(p_values))){
    warning('NA values detected in p-values')
    p_values <- p_values[!is.na(p_values)]
  }
  
  if (length(p_values) == 0) {
    warning('No valid p-values detected')
    return(NA)  # Return NA if no valid p-values
  }
  
  if (method %in% c('bf', "bonferroni",'BF','Bonferroni')) {
    combined_p <- min(1, length(p_values) * min(p_values))
    
  } else if (method %in% c('Fisher',"fisher")) { # Compute Fisher’s combined test statistic
    X2 <- -2 * sum(log(p_values))
    df <- 2 * length(p_values)  
    combined_p <- pchisq(X2, df = df, lower.tail = FALSE)  
  }
  return(combined_p)
}

discretize_fd <- function(x) {
  bin_width <- 2 * IQR(x, na.rm = TRUE) / (length(x)^(1/3))  # Freedman-Diaconis bin width
  nbins <- max(2, round((max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) / bin_width))  # Ensure at least 2 bins
  return(discretize(x, nbins = nbins))
}


# ----------------------------
#  Define the Simulation Function
# ----------------------------

# Function to perform a single simulation iteration
simulate_iteration <- function(cond, iter){
  
  # Extract condition parameters
  #cond<-list(n=500, nvars=4, misfit_type='none') ##run from this line to diagnose simulation
  n <- cond$n
  misfit_type <- cond$misfit_type
  nvars <- cond$nvars
  out <- data.table()
  
  # -----------------------------------
  # Data Generation
  # -----------------------------------
  
  data <- data.table(mvrnorm(n, mu = rep(.5, nvars), Sigma = diag(2,nvars)))
  
  
  data[,V3:=  0 * V1 + 
      0.3 * V2 + 
      ifelse('interaction' %in% misfit_type, .5, 0) * (V1 * V2) + # Introduce nonlinearity (if needed) by an interaction term
      rnorm(n, mean = 0, sd = 1)
  ]
  
  if('groupvar' %in% misfit_type) data[1:floor(n/2),c('V2','V3') ] <- data[1:(n/2),c('V2','V3'),with=F ] * 2
  if('groupmean' %in% misfit_type) data[1:floor(n/2),c('V2','V3') ] <- data[1:(n/2),c('V2','V3'),with=F ] + 2
  
  if('nonnormal' %in% misfit_type) data[,V4:=log1p(exp(V4))]
  
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
      free = ifelse('classical' %in% misfit_type,FALSE,TRUE), values = 0.05, labels = "b2"),
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
  out$SEM_p <-  mxCompare(refs[[1]], fit)$p[2]
  
  # -----------------------------------
  #  Extract Residuals and Perform Misfit Tests
  # -----------------------------------
  
  # Extract residuals
  expMean = c(attributes(fit$output$algebras$PathModel.fitfunction)$expMean)
  expCov <- attributes(fit$output$algebras$PathModel.fitfunction)$expCov
  residuals <- t(apply(data, 1, function(row) row - expMean))
  
  # # Extract row-wise sum of squared standardized residuals
  # stdresidualsrowsumsq <- attributes(fit$output$algebras$PathModel.fitfunction)$rowDist
  
  # compute conditional predictions and residuals
    stdcondresiduals <- matrix(NA, nrow = nrow(data), ncol = nvars)
    condpredictions <- condresiduals <- matrix(NA, nrow = nrow(data), ncol = nvars)
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
      condpredictions[,j] <- expMean[j] + as.vector(Sigma_j_rest %*% Sigma_inv %*% t(X_rest - expMean[-j]))  # Convert to vector

      # Compute residuals
      condresiduals[, j] <- (data[[ colnames(data)[j]]] - condpredictions[,j])
      stdcondresiduals[, j] <- condresiduals[, j] / c(sqrt(cond_var))
    }

  
  stdresiduals <- {
    # 1) Center the columns by the implied means
    # 2) Apply the inverse Cholesky factor to "whiten" them
    t(
      solve(
        chol(expCov),
        t( sweep(as.matrix(data), 2, expMean, FUN = "-") )
      )
    )
  }
  
  # print(cor(stdresiduals))
  
  stdresidualsrowsumsq <- apply(stdresiduals,1,function(x) sum((x^2)))
  stdresidualsrowsumsqsigned <- apply(stdresiduals,1,function(x) sum((x^2)*sign(x)))
  
  
  ##manually compute squared standardized residuals rowsum using expectations
  # stdresiduals <- apply(residuals,1,function(x) x %*% solve(expCov) %*% x)
  
  df <- length(manifests)  # Degrees of freedom for ks chi-square test (number of observed variables)
  
  # Perform tests of std residuals against chi-square distribution
  out$AD_p <- ad.test(stdresidualsrowsumsq, 'pchisq', df = df)$p.value
  # out$ADcondsigned_p <- ad.test(stdresidualsrowsumsqsigned, null = function(x) ecdf(generate_signed_chisq(100000,df))(x))$p.value

  # out$ADboot_p <- {
  #   adfunc <- function(x) ad.test(apply(x,1,function(xx) sum(xx^2)), 'pchisq', df = df)$statistic
  #   observed_ad_stat <- adfunc(stdresiduals) # Compute observed AD test statistic
  #   # Generate null distribution by shuffling within rows
  #   null_ad_stats <- replicate(nboot, {
  #     shuffled_residuals <- matrix(sample(c(stdresiduals)),n,nvars)  # Shuffle within each row
  #     adfunc(shuffled_residuals)  # Compute new ad statistic on shuffled residuals
  #   })
  #   # Compute empirical p-value
  #   mean(null_ad_stats >= observed_ad_stat)
  # }

  
  
  
  
  # 
  # out$CVM_p <-cvm.test(stdresidualsrowsumsq, 'pchisq', df = df)$p.value
  # 
  # out$KS_p <- ks.test(stdresidualsrowsumsq, 'pchisq', df = df)$p.value
  # 
  out$ADuni <- {
    adtestp <- c()
    for(vari in 1:nvars){
      adtestp <- c(adtestp, ad.test(stdresiduals[,vari], 'pnorm',0,1)$p.value)
    }
    combine_pvalues(adtestp) #fisher/ bonferroni correction
  }
  # 
  # out$univarADcondStacked_p <- {
  #   adtest <- c()
  #   stackedstdcondresiduals<-c()
  #   for(vari in 1:nvars){
  #     stackedstdcondresiduals <- c(
  #       stackedstdcondresiduals,
  #       stdresiduals[, vari]
  #     )
  #   }
  #   ad.test(stackedstdcondresiduals, 'pnorm',0,1)$p.value
  # }
  
  
  # out$MVnorm_p <- mvnorm.etest(residuals,nboot)$p.value
  out$MVnormE_p <- mvnorm.etest(stdresiduals,nboot)$p.value

  out$dcov_p <- {   # Test pairwise independence among whitened residuals
    pout <- c()
    for (i in 1:(nvars-1)) {
      for (j in (i+1):nvars) {
        pout <- c(pout,
          dcov.test(stdresiduals[, i],
            stdresiduals[, j],R=nboot)$p.value)
      }
    }
    combine_pvalues(pout)
  }
  
  # #Julian
  # out$dcovpred_p <- {   # Test pairwise independence among predictions and errors
  #   pout <- c()
  #   for (i in 1:nvars) {
  #     for (j in 1:nvars) {
  #       pout <- c(pout,
  #         dcov.test(condpredictions[, i],
  #           condresiduals[, j],R=nboot)$p.value)
  #     }
  #   }
  #   combine_pvalues(pout)
  # }
  
  out$dHSIC_p <- {   # Test pairwise independence among columns
    # pout <- c()
    # for (i in 1:(nvars-1)) {
    #   for (j in (i+1):nvars) {
    #     pout <- c(pout,
    #       dhsic.test(stdresiduals[, i], 
    #         stdresiduals[, j])$p.value)
    #   }
    # }
    dhsic.test(stdresiduals,matrix.input=T,method='gamma')$p.value
    # combine_pvalues(pout)
  }
  
  out$dHSICbiv_p <- {   # Test pairwise independence among columns
    pout <- c()
    for (i in 1:(nvars-1)) {
      for (j in (i+1):nvars) {
        pout <- c(pout,
          dhsic.test(stdresiduals[, i],
            stdresiduals[, j],method='gamma')$p.value)
      }
    }
    combine_pvalues(pout)
  }
  
  
  # # Levene’s test for heteroscedasticity in residual variance
  # Levene_p <- {
  #   levenetest <- c()
  #   for(vari in 1:nvars){
  #     # Grouping continuous predictors into quantiles
  #     group_var <- cut(residuals[, vari], breaks = 3)  
  #     levenetest <- c(levenetest, leveneTest(stdresidualsrowsumsq ~ group_var)$"Pr(>F)"[1])
  #   }
  #   min(levenetest) * nvars # Bonferroni correction
  # }
  
  
  
  
  
  # Mutual Information test for nonlinear dependencies
  out$MIagg_p <- {
    mi_results <- c()
    p_values <- c()
    for (vari in 1:nvars) {
      stdcondresidualsrowsumrest <- apply(stdresiduals[, -vari,drop=F],1,function(x) sum((x)))
      discCondResiduals <- discretize(stdresiduals[, vari])
      # Compute observed MI
      observed_mi <- mutinformation(discretize(stdcondresidualsrowsumrest), discCondResiduals)
      # Generate null distribution by permuting R
      permuted_mis <- replicate(nboot, {
        mutinformation(discretize(sample(stdcondresidualsrowsumrest)), discCondResiduals)
      })
      # Compute p-value: Proportion of permuted MIs greater than or equal to observed MI
      p_values <- c(p_values, mean(permuted_mis >= observed_mi))
    }
    # Apply correction for multiple comparisons
    combine_pvalues(p_values)
  }
  
  # out$MutInfoCondRest2_p <- {
  #   mi_results <- c()
  #   p_values <- c()
  #   for (vari in 1:nvars) {
  #     stdcondresidualsrowsumrest <- apply(stdresiduals[, -vari,drop=F],1,function(x) sum((x^2)*sign(x)))
  #     discCondResiduals <- discretize(stdresiduals[, vari]^2*sign(stdresiduals[, vari]))
  #     # Compute observed MI
  #     observed_mi <- mutinformation(discretize(stdcondresidualsrowsumrest), discCondResiduals)
  #     # Generate null distribution by permuting R
  #     permuted_mis <- replicate(nboot, {
  #       mutinformation(discretize(sample(stdcondresidualsrowsumrest)), discCondResiduals)
  #     })
  #     # Compute p-value: Proportion of permuted MIs greater than or equal to observed MI
  #     p_values <- c(p_values, mean(permuted_mis >= observed_mi))
  #   }
  #   # Apply correction for multiple comparisons
  #   combine_pvalues(p_values)
  # }
  
  # out$MutInfoCondRest2abs_p <- {
  #   mi_results <- c()
  #   p_values <- c()
  #   for (vari in 1:nvars) {
  #     stdcondresidualsrowsumrest <- apply(stdresiduals[, -vari,drop=F],1,function(x) sum((x^2)))
  #     discCondResiduals <- discretize(stdresiduals[, vari]^2)
  #     # Compute observed MI
  #     observed_mi <- mutinformation(discretize(stdcondresidualsrowsumrest), discCondResiduals)
  #     # Generate null distribution by permuting R
  #     permuted_mis <- replicate(nboot, {
  #       mutinformation(discretize(sample(stdcondresidualsrowsumrest)), discCondResiduals)
  #     })
  #     # Compute p-value: Proportion of permuted MIs greater than or equal to observed MI
  #     p_values <- c(p_values, mean(permuted_mis >= observed_mi))
  #   }
  #   # Apply correction for multiple comparisons
  #   combine_pvalues(p_values)
  # }
  
  # out$MutInfoCondRest2absShuffle_p <- {
  #  mi_results <- c()
  #  p_values <- c()
  #  for (vari in 1:nvars) {
  #    stdcondresidualsrowsumrest <- apply(stdresiduals[, -vari,drop=F],1,function(x) sum((x^2)))
  #    discCondResiduals <- discretize(stdresiduals[, vari])
  #    # Compute observed MI
  #    observed_mi <- mutinformation(discretize(stdcondresidualsrowsumrest), discCondResiduals)
  #    # Generate null distribution by permuting R
  #    permuted_mis <- replicate(nboot, {
  #      mutinformation(
  #        discretize(sample(apply(matrix(sample(c(stdresiduals[, -vari]),replace=FALSE),n,nvars-1),1,function(x) sum((x^2))))), 
  #        discCondResiduals)
  #    })
  #    # Compute p-value: Proportion of permuted MIs greater than or equal to observed MI
  #    p_values <- c(p_values, mean(permuted_mis >= observed_mi))
  #  }
  #  # Apply correction for multiple comparisons
  #  combine_pvalues(p_values)
  # }
  
  # out$MutInfoCondRest2absShuffle2_p <- {
  #   mi_results <- c()
  #   p_values <- c()
  #   for (vari in 1:nvars) {
  #     stdcondresidualsrowsumrest <- apply(stdresiduals[, -vari,drop=F],1,function(x) sum((x^2)))
  #     discCondResiduals <- discretize(stdresiduals[, vari]^2)
  #     # Compute observed MI
  #     observed_mi <- mutinformation(discretize(stdcondresidualsrowsumrest), discCondResiduals)
  #     # Generate null distribution by permuting R
  #     permuted_mis <- replicate(nboot, {
  #       mutinformation(
  #         discretize(sample(apply(matrix(sample(c(stdresiduals[, -vari]),replace=FALSE),n,nvars-1),1,function(x) sum((x^2))))), 
  #         discCondResiduals)
  #     })
  #     # Compute p-value: Proportion of permuted MIs greater than or equal to observed MI
  #     p_values <- c(p_values, mean(permuted_mis >= observed_mi))
  #   }
  #   # Apply correction for multiple comparisons
  #   combine_pvalues(p_values)
  # }
  
  out$MIbiv_p <- {
    mi_results <- c()
    p_values <- c()

    for (var1 in 1:(nvars-1)) {  # Iterate over first variable in pair
      for (var2 in (var1+1):nvars) {  # Iterate over second variable in pair

        # Extract conditional residuals for the two variables
        discRes1 <- discretize(stdresiduals[, var1])
        discRes2 <- discretize(stdresiduals[, var2])

        # Compute observed MI
        observed_mi <- mutinformation(discRes1, discRes2)

        # Generate null distribution by permuting one variable independently
        permuted_mis <- replicate(nboot, {
          shuffled_res <- discretize(sample(stdresiduals[, var2]))
          mutinformation(discRes1, shuffled_res)
        })

        # Compute p-value: Proportion of permuted MIs greater than or equal to observed MI
        p_values <- c(p_values, mean(permuted_mis >= observed_mi))
      }
    }
    combine_pvalues(p_values)
  }
  
  
  # p <- c()
  # for(i in 1:100){
  # stdresiduals <- mvrnorm(n, mu = rep(0, nvars), Sigma = diag(1,nvars))
  # 
  # out$MutInfoCondRestStacked_p <- {
  #   mi_results <- c()
  #   p_values <- c()
  #   stackedstdcondresidualsrowsumrest<-c()
  #   stackedCondResiduals<-c()
  #   for (vari in 1:nvars) {
  #     stackedstdcondresidualsrowsumrest <- c(
  #       stackedstdcondresidualsrowsumrest,
  #       apply(stdresiduals[, -vari,drop=F],1,function(x) sum((x^2)*sign(x)))
  #     )
  #     stackedCondResiduals <- c(
  #       stackedCondResiduals,
  #       stdresiduals[, vari]
  #     )
  #   }
  #   nbins <- (n)^(1/3)
  #   discStackedCondResiduals <- discretize(stackedCondResiduals,nbins=nbins)
  #   # Compute observed MI
  #   observed_mi <- mutinformation(discretize(stackedstdcondresidualsrowsumrest,nbins=nbins), discStackedCondResiduals)
  #   # Generate null distribution by permuting R
  #   permuted_mis <- replicate(nboot, {
  #     mutinformation(discretize(sample(stackedstdcondresidualsrowsumrest),nbins=nbins), 
  #     discStackedCondResiduals)
  #   })
  #   # Compute p-value: Proportion of permuted MIs greater than or equal to observed MI
  #   p_value<- mean(permuted_mis >= observed_mi)
  # }
  # p <- c(p,p_value)
  # }
  # print(mean(p<.05))
  
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
      a <- 0.032  # Scaling factorp
      return((1 / (a * h_het + 1))-.9)  # HFI output less than .05 = significant, -.9 was added by me. 
    }
  }
  
  # out$HFIunstd=compute_HFI(residuals)
  out$HFI=compute_HFI(stdresiduals)
  
  # Return p-values and HFI
  return(out)
}

# ----------------------------
#  Run Simulations in Parallel
# ----------------------------

# Use foreach to iterate over each simulation condition
# Outer loop iterates over conditions
# Inner loop iterates over iterations within each condition
simresultslist <- foreach(cond_idx = 1:nrow(simconditions), #.combine = rbind, 
  .packages = c("OpenMx", "MASS", "data.table", "moments", "goftest",'energy', 'infotheo','dHSIC')) %:%
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
outputstatnames <- colnames(simresults)[!colnames(simresults) %in% c('condition','n','misfit_type','nvars')]
summarynames <- gsub('_p','',outputstatnames) #remove _p from names
summary_results <- simresults[, lapply(.SD, function(x) round(mean(x <.05, na.rm = TRUE), 3)), by = .(condition, n, nvars, misfit_type), .SDcols = outputstatnames]

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


summary_long[,Misfit:=misfit_type]
summary_long[,misfit_metric:=gsub('_p','',misfit_metric)]
summary_long[,n:=factor(n)]
# Plot misfit rates for SEM chi-square test, Raw KS test, and HFI across conditions
p <- ggplot(
  summary_long[grepl('(MV)|(HFI)|(HSIC)|(SEM)|(dcov)|(AD)|(MI)',misfit_metric),],#[!misfit_metric %in% c('KS','CVM','Levene'),],
  aes(x = misfit_metric, y = misfit_rate, fill = misfit_metric)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  facet_grid(rows = vars(Misfit),cols=vars(interaction(nvars,n)),labeller = label_both) +
  labs(
    # x = "Sample Size (n)",
    y = "Proportion of p-values ≤ 0.05",
    fill = "Misfit Metric"
  ) +
  theme_bw() +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))

print(p)
if(T){
  pdf('sim.pdf', width = 10, height = 12)
  print(p)
  dev.off()
}
# ================================================================
# End of Simulation Study Script
# ================================================================
