# Full WP10B run; checkpoints let `Rscript WP10B_calibration_simulation.R` resume safely.
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
library(data.table); library(MASS); library(Matrix); library(ggplot2); library(kableExtra)
source("sem_misfit_api_prototype.R")
seed_master <- 20260721L
set.seed(seed_master)
B <- 499L; ncores <- 10L; alpha <- .05
conditions <- c("gaussian_null", "marginal_skewness_only", "heavy_tail_only", "covariance_misspecification", "nonlinear_conditional_mean", "conditional_heteroscedasticity", "subgroup_heterogeneity", "rotated_independent_nongaussian", "higher_order_only_dependence")
comparison_conditions <- c("gaussian_null", "covariance_misspecification", "nonlinear_conditional_mean", "conditional_heteroscedasticity")
simulation_spec <- "wp10b_refit_primary_B499_familywise_v3"
output_dir <- paste0("wp10b_results_", simulation_spec)
if (identical(Sys.getenv("WP10B_SMOKE"), "1")) output_dir <- paste0(output_dir, "_smoke")
checkpoint_dir <- file.path(output_dir, "checkpoints")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
design <- CJ(n=c(100L,500L), nvars=c(4L,8L), condition=conditions)
design[, outer_reps := ifelse(condition == "gaussian_null", 1000L, 500L)]
design[, cell_id := sprintf("n%d_p%d_%s",n,nvars,condition)]
if (identical(Sys.getenv("WP10B_SMOKE"), "1")) { design <- design[1L]; design[, outer_reps := 1L] }

fit_model <- function(data) {
 x <- as.matrix(data); p <- ncol(x); s <- cov(x)*(nrow(x)-1)/nrow(x); sigma <- diag(diag(s),p)
 for (ij in list(c(1L,3L),c(2L,3L))) sigma[ij[1],ij[2]] <- sigma[ij[2],ij[1]] <- s[ij[1],ij[2]]
 sigma <- as.matrix(nearPD(sigma,keepDiag=TRUE)$mat)
 as_sem_expectation(list(data=as.data.frame(data),mean=colMeans(x),covariance=sigma,variables=names(data),source="wp10b_moment_model",refit=function(expectation,data) fit_model(data)))
}
make_data <- function(n,p,condition) {
 x <- mvrnorm(n,rep(0,p),diag(p)); colnames(x) <- paste0("V",seq_len(p)); x[,3] <- .35*x[,1]+.4*x[,2]+rnorm(n)
 if(condition=="marginal_skewness_only") x[,4] <- scale(exp(x[,4]))[,1]
 if(condition=="heavy_tail_only") x[,4] <- rt(n,3)/sqrt(3)
 if(condition=="covariance_misspecification") x[,1] <- x[,1]+.55*x[,4]
 if(condition %in% c("nonlinear_conditional_mean","higher_order_only_dependence")) x[,3] <- x[,3]+.65*(x[,1]^2-1)
 if(condition=="conditional_heteroscedasticity") x[,3] <- .35*x[,1]+.4*x[,2]+(1+.8*abs(x[,1]))*rnorm(n)
 if(condition=="subgroup_heterogeneity") {g<-rbinom(n,1,.5);x[,3]<-x[,3]+g*(1.1+.6*x[,1])}
 if(condition=="rotated_independent_nongaussian") {u<-apply(matrix(rexp(n*p)-1,n,p),2,scale);x<-u%*%qr.Q(qr(matrix(rnorm(p*p),p)));colnames(x)<-paste0("V",seq_len(p))}
 as.data.frame(x)
}
targets <- function(condition) {
 if(condition %in% c("nonlinear_conditional_mean", "higher_order_only_dependence")) return(list(variable="V3", pairs=c("V3|V1", "V3|V2"), scale=NA_character_))
 if(condition == "conditional_heteroscedasticity") return(list(variable=NA_character_, pairs=character(), scale="V3"))
 if(condition == "covariance_misspecification") return(list(variable=NA_character_, pairs=c("V1|V4", "V4|V1"), scale=NA_character_))
 list(variable=NA_character_, pairs=character(), scale=NA_character_)
}
calibration_modes <- function(condition) {
  if (condition %in% comparison_conditions) c("conditional_permutation", "parametric_fixed", "parametric_refit") else "parametric_refit"
}
format_elapsed <- function(seconds) {
 seconds <- max(0, round(seconds)); sprintf("%02d:%02d:%02d", seconds %/% 3600L, (seconds %% 3600L) %/% 60L, seconds %% 60L)
}
family_reject <- function(tab) any(!is.na(tab$p_adjusted) & tab$p_adjusted <= alpha)
run_one <- function(row,i) {
 set.seed(seed_master+as.integer(i+10000*match(row$cell_id,design$cell_id))); e<-fit_model(make_data(row$n,row$nvars,row$condition)); target<-targets(row$condition)
 out <- rbindlist(lapply(calibration_modes(row$condition),function(mode){
  started<-proc.time()[["elapsed"]]; z<-tryCatch(sem_misfit(e,tests=c("conditional_variables","conditional_pairs","conditional_scale"),calibration=mode,B=B,seed=sample.int(.Machine$integer.max,1)),error=function(err)err)
  if(inherits(z,"error")) return(data.table(calibration=mode,variable_family_reject=NA,pair_family_reject=NA,scale_family_reject=NA,any_family_signal=NA,variable_correct=NA,pair_correct=NA,scale_correct=NA,runtime=proc.time()[["elapsed"]]-started,success=0L,failure_rate=1,error=conditionMessage(z)))
  cv<-z$conditional_variables; cp<-z$conditional_pairs; cs<-z$conditional_scale; keys<-paste(cp$residual_variable,cp$predictor,sep="|")
  variable_reject <- family_reject(cv); pair_reject <- family_reject(cp); scale_reject <- family_reject(cs)
  data.table(calibration=mode,variable_family_reject=variable_reject,pair_family_reject=pair_reject,scale_family_reject=scale_reject,any_family_signal=any(variable_reject,pair_reject,scale_reject),variable_correct=if(is.na(target$variable)) NA else any(cv$variable==target$variable & !is.na(cv$p_adjusted) & cv$p_adjusted<=alpha),pair_correct=if(!length(target$pairs)) NA else any(keys%in%target$pairs & !is.na(cp$p_adjusted) & cp$p_adjusted<=alpha),scale_correct=if(is.na(target$scale)) NA else any(cs$variable==target$scale & !is.na(cs$p_adjusted) & cs$p_adjusted<=alpha),runtime=z$calibration_details$runtime_seconds%||%(proc.time()[["elapsed"]]-started),success=z$calibration_details$B_successful%||%NA_integer_,failure_rate=if(is.null(z$calibration_details))NA_real_ else 1-z$calibration_details$B_successful/z$calibration_details$B_requested,error=NA_character_)
 }))
 out[, `:=`(n = row$n, nvars = row$nvars, condition = row$condition, iter = i, simulation_spec = simulation_spec)]
 out
}
run_cell <- function(row) {
 path <- file.path(checkpoint_dir, paste0(row$cell_id, ".rds"))
 if (file.exists(path)) {
   cached <- tryCatch(readRDS(path), error = function(e) e)
   expected_rows <- length(calibration_modes(row$condition)) * row$outer_reps
   if (!inherits(cached, "error") && nrow(cached) == expected_rows && "simulation_spec" %in% names(cached) && all(cached$simulation_spec == simulation_spec)) return(cached)
   message("Ignoring incomplete or unreadable checkpoint for ", row$cell_id, "; recomputing this cell.")
 }
 message("Running ",row$cell_id,"; outer replications: ",row$outer_reps)
 cl<-parallel::makeCluster(ncores);on.exit(parallel::stopCluster(cl),add=TRUE)
 parallel::clusterEvalQ(cl,{.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()));library(data.table);library(MASS);library(Matrix);library(energy);library(goftest);library(dHSIC);source("sem_misfit_api_prototype.R")})
 parallel::clusterExport(cl,c("seed_master","B","alpha","design","conditions","comparison_conditions","simulation_spec","fit_model","make_data","targets","calibration_modes","family_reject","run_one","%||%"),envir=environment());parallel::clusterSetRNGStream(cl,seed_master+match(row$cell_id,design$cell_id))
 progress_every <- max(ncores, ceiling(row$outer_reps / 20L))
 batches <- split(seq_len(row$outer_reps), ceiling(seq_len(row$outer_reps) / progress_every))
 cell_started <- proc.time()[["elapsed"]]; completed <- 0L; pieces <- vector("list", length(batches))
 for (b in seq_along(batches)) {
   pieces[[b]] <- parallel::parLapplyLB(cl, batches[[b]], function(i) run_one(row, i))
   completed <- completed + length(batches[[b]]); elapsed <- proc.time()[["elapsed"]] - cell_started
   remaining <- if (completed) elapsed / completed * (row$outer_reps - completed) else NA_real_
   message(sprintf("%s: %d/%d (%.1f%%), elapsed %s, estimated remaining %s", row$cell_id, completed, row$outer_reps, 100 * completed / row$outer_reps, format_elapsed(elapsed), format_elapsed(remaining))); flush.console()
 }
 ans <- rbindlist(unlist(pieces, recursive = FALSE), fill = TRUE)
 temporary_path <- tempfile(pattern = paste0(row$cell_id, "-"), tmpdir = checkpoint_dir, fileext = ".rds")
 saveRDS(ans, temporary_path)
 if (!file.copy(temporary_path, path, overwrite = TRUE)) stop("Could not write checkpoint: ", path, call. = FALSE)
 unlink(temporary_path)
 ans
}
started<-Sys.time();simresults<-rbindlist(lapply(seq_len(nrow(design)),function(i)run_cell(design[i])),fill=TRUE);save(simresults,design,B,ncores,seed_master,simulation_spec,file=file.path(output_dir,"simresults.rda"))
wilson<-function(k,n){if(!n)return(c(NA_real_,NA_real_));z<-qnorm(.975);p<-k/n;d<-1+z^2/n;c((p+z^2/(2*n)-z*sqrt(p*(1-p)/n+z^2/(4*n^2)))/d,(p+z^2/(2*n)+z*sqrt(p*(1-p)/n+z^2/(4*n^2)))/d)}
family_columns <- c(variable_family_reject="variable_family",pair_family_reject="directional_pair_family",scale_family_reject="conditional_scale_family",any_family_signal="any_family_signal_descriptive")
summary<-rbindlist(lapply(names(family_columns),function(column)simresults[,.(diagnostic_family=family_columns[[column]],rejection_rate=mean(get(column),na.rm=TRUE),mc_low=wilson(sum(get(column),na.rm=TRUE),sum(!is.na(get(column))))[1],mc_high=wilson(sum(get(column),na.rm=TRUE),sum(!is.na(get(column))))[2],median_runtime=median(runtime,na.rm=TRUE),bootstrap_failure_rate=mean(failure_rate,na.rm=TRUE),localization_accuracy=mean(get(switch(column,variable_family_reject="variable_correct",pair_family_reject="pair_correct",scale_family_reject="scale_correct",any_family_signal="variable_correct")),na.rm=TRUE),unavailable=sum(!is.na(error))),by=.(n,nvars,condition,calibration)]),fill=TRUE)
fwrite(summary,file.path(output_dir,"wp10b_calibration_summary.csv"));sink(file.path(output_dir,"simtable.tex"));print(kableExtra::kable(as.data.frame(summary),format="latex",booktabs=TRUE,digits=3));sink()
p<-ggplot(summary[diagnostic_family!="any_family_signal_descriptive"],aes(condition,rejection_rate,fill=calibration))+geom_col(position="dodge")+geom_errorbar(aes(ymin=mc_low,ymax=mc_high),position=position_dodge(.9),width=.2)+facet_grid(diagnostic_family+nvars~n,labeller=label_both)+coord_cartesian(ylim=c(0,1))+theme_bw(base_size=8)+theme(axis.text.x=element_text(angle=35,hjust=1))+labs(x=NULL,y="Family-wise rejection rate (95% binomial MC interval)",fill="Calibration")
ggsave(file.path(output_dir,"sim.pdf"),p,width=13,height=12);packages<-c("data.table","MASS","Matrix","ggplot2","kableExtra","energy","goftest","dHSIC");writeLines(c(paste("started",started),paste("finished",Sys.time()),paste("simulation_spec",simulation_spec),paste("B",B),paste("cores",ncores),paste("seed_master",seed_master),paste("comparison conditions",paste(comparison_conditions,collapse=", ")),paste("package versions",paste(paste0(packages,"=",vapply(packages,function(x)as.character(packageVersion(x)),character(1))),collapse="; ")),"family rejection uses p_adjusted; any-family signal is descriptive only",paste("checkpoint directory",checkpoint_dir),paste("output directory",output_dir)),file.path(output_dir,"wp10b_run_metadata.txt"))
