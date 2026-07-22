# Full WP10B run; checkpoints let `Rscript WP10B_calibration_simulation.R` resume safely.
library(data.table); library(MASS); library(Matrix); library(ggplot2); library(kableExtra)
source("sem_misfit_api_prototype.R")
set.seed(20260721)
B <- 999L; ncores <- 10L; alpha <- .05; checkpoint_dir <- "wp10b_checkpoints"
dir.create(checkpoint_dir, showWarnings = FALSE)
conditions <- c("gaussian_null", "marginal_skewness_only", "heavy_tail_only", "covariance_misspecification", "nonlinear_conditional_mean", "conditional_heteroscedasticity", "subgroup_heterogeneity", "rotated_independent_nongaussian", "higher_order_only_dependence")
controls <- conditions[1:3]
design <- CJ(n=c(100L,500L), nvars=c(4L,8L), condition=conditions)
design[, outer_reps := ifelse(condition %in% controls, 1000L, 500L)]
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
targets <- function(condition) if(condition %in% c("nonlinear_conditional_mean","conditional_heteroscedasticity","higher_order_only_dependence")) list(variable="V3",pairs=c("V3|V1","V3|V2")) else if(condition=="covariance_misspecification") list(variable=NA_character_,pairs=c("V1|V4","V4|V1")) else list(variable=NA_character_,pairs=character())
format_elapsed <- function(seconds) {
 seconds <- max(0, round(seconds)); sprintf("%02d:%02d:%02d", seconds %/% 3600L, (seconds %% 3600L) %/% 60L, seconds %% 60L)
}
run_one <- function(row,i) {
 set.seed(20260721L+as.integer(i+10000*match(row$cell_id,design$cell_id))); e<-fit_model(make_data(row$n,row$nvars,row$condition)); target<-targets(row$condition)
 rbindlist(lapply(c("conditional_permutation","parametric_fixed","parametric_refit"),function(mode){
  started<-proc.time()[["elapsed"]]; z<-tryCatch(sem_misfit(e,tests=c("conditional_variables","conditional_pairs","conditional_scale"),calibration=mode,B=B,seed=sample.int(.Machine$integer.max,1)),error=function(err)err)
  if(inherits(z,"error")) return(data.table(calibration=mode,rejection=NA,variable_correct=NA,pair_correct=NA,runtime=proc.time()[["elapsed"]]-started,success=0L,failure_rate=1,error=conditionMessage(z)))
  cv<-z$conditional_variables;cp<-z$conditional_pairs; keys<-paste(cp$residual_variable,cp$predictor,sep="|")
  data.table(calibration=mode,rejection=any(c(cv$p_value,cp$p_value,z$conditional_scale$p_value)<=alpha,na.rm=TRUE),variable_correct=if(is.na(target$variable)) NA else any(cv$variable==target$variable&cv$p_value<=alpha),pair_correct=if(!length(target$pairs)) NA else any(keys%in%target$pairs&cp$p_value<=alpha),runtime=z$calibration_details$runtime_seconds%||%(proc.time()[["elapsed"]]-started),success=z$calibration_details$B_successful%||%NA_integer_,failure_rate=if(is.null(z$calibration_details))NA_real_ else 1-z$calibration_details$B_successful/z$calibration_details$B_requested,error=NA_character_)
 }))[,`:=`(n=row$n,nvars=row$nvars,condition=row$condition,iter=i)]
}
run_cell <- function(row) {
 path<-file.path(checkpoint_dir,paste0(row$cell_id,".rds")); if(file.exists(path)) { cached<-readRDS(path); if(nrow(cached)==3L*row$outer_reps)return(cached) }; message("Running ",row$cell_id,"; outer replications: ",row$outer_reps)
 cl<-parallel::makeCluster(ncores);on.exit(parallel::stopCluster(cl),add=TRUE)
 parallel::clusterEvalQ(cl,{library(data.table);library(MASS);library(Matrix);library(energy);library(goftest);library(dHSIC);source("sem_misfit_api_prototype.R")})
 parallel::clusterExport(cl,c("B","alpha","design","conditions","fit_model","make_data","targets","run_one","%||%"),envir=environment());parallel::clusterSetRNGStream(cl,20260721L+match(row$cell_id,design$cell_id))
 progress_every <- max(ncores, ceiling(row$outer_reps / 20L))
 batches <- split(seq_len(row$outer_reps), ceiling(seq_len(row$outer_reps) / progress_every))
 cell_started <- proc.time()[["elapsed"]]; completed <- 0L; pieces <- vector("list", length(batches))
 for (b in seq_along(batches)) {
   pieces[[b]] <- parallel::parLapplyLB(cl, batches[[b]], function(i) run_one(row, i))
   completed <- completed + length(batches[[b]]); elapsed <- proc.time()[["elapsed"]] - cell_started
   remaining <- if (completed) elapsed / completed * (row$outer_reps - completed) else NA_real_
   message(sprintf("%s: %d/%d (%.1f%%), elapsed %s, estimated remaining %s", row$cell_id, completed, row$outer_reps, 100 * completed / row$outer_reps, format_elapsed(elapsed), format_elapsed(remaining))); flush.console()
 }
 ans<-rbindlist(unlist(pieces,recursive=FALSE),fill=TRUE);saveRDS(ans,path);ans
}
started<-Sys.time();simresults<-rbindlist(lapply(seq_len(nrow(design)),function(i)run_cell(design[i])),fill=TRUE);save(simresults,design,B,ncores,file="simresults.rda")
wilson<-function(k,n){z<-qnorm(.975);p<-k/n;d<-1+z^2/n;c((p+z^2/(2*n)-z*sqrt(p*(1-p)/n+z^2/(4*n^2)))/d,(p+z^2/(2*n)+z*sqrt(p*(1-p)/n+z^2/(4*n^2)))/d)}
summary<-simresults[,.(rejection_rate=mean(rejection,na.rm=TRUE),mc_low=wilson(sum(rejection,na.rm=TRUE),sum(!is.na(rejection)))[1],mc_high=wilson(sum(rejection,na.rm=TRUE),sum(!is.na(rejection)))[2],median_runtime=median(runtime,na.rm=TRUE),bootstrap_failure_rate=mean(failure_rate,na.rm=TRUE),variable_accuracy=mean(variable_correct,na.rm=TRUE),pair_accuracy=mean(pair_correct,na.rm=TRUE),unavailable=sum(!is.na(error))),by=.(n,nvars,condition,calibration)]
fwrite(summary,"wp10b_calibration_summary.csv");sink("simtable.tex");print(kableExtra::kable(as.data.frame(summary),format="latex",booktabs=TRUE,digits=3));sink()
p<-ggplot(summary,aes(condition,rejection_rate,fill=calibration))+geom_col(position="dodge")+geom_errorbar(aes(ymin=mc_low,ymax=mc_high),position=position_dodge(.9),width=.2)+facet_grid(nvars~n,labeller=label_both)+coord_cartesian(ylim=c(0,1))+theme_bw(base_size=9)+theme(axis.text.x=element_text(angle=35,hjust=1))+labs(x=NULL,y="Rejection rate (95% binomial Monte Carlo interval)",fill="Calibration")
ggsave("sim.pdf",p,width=13,height=9);writeLines(c(paste("started",started),paste("finished",Sys.time()),paste("B",B),paste("cores",ncores),"checkpoint directory: wp10b_checkpoints"),"wp10b_run_metadata.txt")
