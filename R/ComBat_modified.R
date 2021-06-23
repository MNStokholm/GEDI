## THIS SCRIPT USES MODIFIED FUNCTIONS FROM THE SVA PACKAGE
## https://bioconductor.org/packages/release/bioc/html/sva.html
## MODIFIED LINES ARE CLEARLY MARKED WITH COMMENTS WITH THE WORD 'MODIFIED'


#### Helper #################################################################
rowVars = function(x, ...) {
  sqr     = function(x)  x*x
  n       = rowSums(!is.na(x))
  n[n<=1] = NA
  return(rowSums(sqr(x-rowMeans(x, ...)), ...)/(n-1))
}

# Following four find empirical hyper-prior values
aprior <- function(gamma.hat) {
  m <- mean(gamma.hat)
  s2 <- var(gamma.hat)
  (2*s2 + m^2) / s2
}

bprior <- function(gamma.hat){
  m <- mean(gamma.hat)
  s2 <- var(gamma.hat)
  (m*s2 + m^3) / s2
}

postmean <- function(g.hat,g.bar,n,d.star,t2){
  (t2*n*g.hat + d.star*g.bar) / (t2*n + d.star)
}

postvar <- function(sum2,n,a,b){
  (.5*sum2 + b) / (n/2 + a - 1)
}

# Inverse gamma distribution density function. (Note: does not do any bounds checking on arguments)
dinvgamma <- function (x, shape, rate = 1/scale, scale = 1) {
  # PDF taken from https://en.wikipedia.org/wiki/Inverse-gamma_distribution
  # Note: alpha = shape, beta = rate
  stopifnot(shape > 0)
  stopifnot(rate > 0)
  ifelse(x <= 0, 0, ((rate ^ shape) / gamma(shape)) * x ^ (-shape - 1) * exp(-rate/x))
}

# Pass in entire data set, the design matrix for the entire data, the batch means, the batch variances, priors (m, t2, a, b), columns of the data  matrix for the batch. Uses the EM to find the parametric batch adjustments

it.sol  <- function(sdat,g.hat,d.hat,g.bar,t2,a,b,conv=.0001){
  n <- rowSums(!is.na(sdat))
  g.old <- g.hat
  d.old <- d.hat
  change <- 1
  count <- 0
  while(change>conv){
    g.new <- postmean(g.hat, g.bar, n, d.old, t2)
    sum2 <- rowSums((sdat - g.new %*% t(rep(1,ncol(sdat))))^2, na.rm=TRUE)
    d.new <- postvar(sum2, n, a, b)
    change <- max(abs(g.new-g.old) / g.old, abs(d.new-d.old) / d.old)
    g.old <- g.new
    d.old <- d.new
    count <- count+1
  }
  ## cat("This batch took", count, "iterations until convergence\n")
  adjust <- rbind(g.new, d.new)
  rownames(adjust) <- c("g.star","d.star")
  adjust
}

## Monte Carlo integration functions
int.eprior <- function(sdat, g.hat, d.hat){
  g.star <- d.star <- NULL
  r <- nrow(sdat)
  for(i in 1:r){
    g <- g.hat[-i]
    d <- d.hat[-i]
    x <- sdat[i,!is.na(sdat[i,])]
    n <- length(x)
    j <- numeric(n)+1
    dat <- matrix(as.numeric(x), length(g), n, byrow=TRUE)
    resid2 <- (dat-g)^2
    sum2 <- resid2 %*% j
    LH <- 1/(2*pi*d)^(n/2)*exp(-sum2/(2*d))
    LH[LH=="NaN"]=0
    g.star <- c(g.star, sum(g*LH)/sum(LH))
    d.star <- c(d.star, sum(d*LH)/sum(LH))
    ## if(i%%1000==0){cat(i,'\n')}
  }
  adjust <- rbind(g.star,d.star)
  rownames(adjust) <- c("g.star","d.star")
  adjust
}

## fits the L/S model in the presence of missing data values

Beta.NA <- function(y,X){
  des <- X[!is.na(y),]
  y1 <- y[!is.na(y)]
  B <- solve(crossprod(des), crossprod(des, y1))
  B
}


#### ComBat ####
#' @importFrom graphics lines par
#' @importFrom stats cor density dnorm model.matrix pf ppoints prcomp predict
#' qgamma qnorm qqline qqnorm qqplot smooth.spline var
#' @importFrom utils read.delim
#' @importFrom BiocParallel bplapply bpparam
.myComBat <- function(dat, batch, mod = NULL, par.prior = TRUE, prior.plots = FALSE,
                   mean.only = FALSE, ref.batch = NULL, BPPARAM = bpparam("SerialParam")) {
  if(length(dim(batch))>1){
    stop("This version of ComBat only allows one batch variable")
  }  ## to be updated soon!

  ## coerce dat into a matrix
  dat <- as.matrix(dat)

  batch <- as.factor(batch)

  ## MODIFIED: Following was removed to correct batch effect in zero variance genes
  # # find genes with zero variance in any of the batches
  # zero.rows.lst <- lapply(levels(batch), function(batch_level){
  #   if(sum(batch==batch_level)>1){
  #     return(which(apply(dat[, batch==batch_level], 1, function(x){var(x)==0})))
  #   }else{
  #     return(which(rep(1,3)==2))
  #   }
  # })
  # zero.rows <- Reduce(union, zero.rows.lst)
  # keep.rows <- setdiff(1:nrow(dat), zero.rows)
  #
  # if (length(zero.rows) > 0) {
  #   cat(sprintf("Found %d genes with uniform expression within a single batch (all zeros); these will not be adjusted for batch.\n", length(zero.rows)))
  #   # keep a copy of the original data matrix and remove zero var rows
  #   dat.orig <- dat
  #   dat <- dat[keep.rows, ]
  # }

  ## make batch a factor and make a set of indicators for batch
  if(any(table(batch)==1)){mean.only=TRUE}
  if(mean.only==TRUE){
    message("Using the 'mean only' version of ComBat")
  }

  batchmod <- model.matrix(~-1+batch)
  if (!is.null(ref.batch)){
    ## check for reference batch, check value, and make appropriate changes
    if (!(ref.batch%in%levels(batch))) {
      stop("reference level ref.batch is not one of the levels of the batch variable")
    }
    message("Using batch =",ref.batch, "as a reference batch (this batch won't change)")
    ref <- which(levels(as.factor(batch))==ref.batch) # find the reference
    batchmod[,ref] <- 1
  } else {
    ref <- NULL
  }
  message("Found", nlevels(batch), "batches")

  ## A few other characteristics on the batches
  n.batch <- nlevels(batch)
  batches <- list()
  for (i in 1:n.batch) {
    batches[[i]] <- which(batch == levels(batch)[i])
  } # list of samples in each batch
  n.batches <- sapply(batches, length)
  if(any(n.batches==1)){
    mean.only=TRUE
    message("Note: one batch has only one sample, setting mean.only=TRUE")
  }
  n.array <- sum(n.batches)
  ## combine batch variable and covariates
  design <- cbind(batchmod,mod)

  ## check for intercept in covariates, and drop if present
  check <- apply(design, 2, function(x) all(x == 1))
  if(!is.null(ref)){
    check[ref] <- FALSE
  } ## except don't throw away the reference batch indicator
  design <- as.matrix(design[,!check])

  ## Number of covariates or covariate levels
  message("Adjusting for", ncol(design)-ncol(batchmod), 'covariate(s) or covariate level(s)')

  ## Check if the design is confounded
  if(qr(design)$rank < ncol(design)) {
    ## if(ncol(design)<=(n.batch)){stop("Batch variables are redundant! Remove one or more of the batch variables so they are no longer confounded")}
    if(ncol(design)==(n.batch+1)) {
      stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat")
    }
    if(ncol(design)>(n.batch+1)) {
      if((qr(design[,-c(1:n.batch)])$rank<ncol(design[,-c(1:n.batch)]))){
        stop('The covariates are confounded! Please remove one or more of the covariates so the design is not confounded')
      } else {
        stop("At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat")
      }
    }
  }

  ## Check for missing values
  NAs <- any(is.na(dat))
  if(NAs){
    message(c('Found',sum(is.na(dat)),'Missing Data Values'), sep=' ')}
  ## print(dat[1:2,])

  ##Standardize Data across genes
  message('Standardizing Data across genes')
  if (!NAs){
    B.hat <- solve(crossprod(design), tcrossprod(t(design), as.matrix(dat)))
  } else {
    B.hat <- apply(dat, 1, Beta.NA, design) # FIXME
  }

  ## change grand.mean for ref batch
  if(!is.null(ref.batch)){
    grand.mean <- t(B.hat[ref, ])
  } else {
    grand.mean <- crossprod(n.batches/n.array, B.hat[1:n.batch,])
  }

  ## change var.pooled for ref batch
  if (!NAs){
    if(!is.null(ref.batch)) {
      ref.dat <- dat[, batches[[ref]]]
      var.pooled <- ((ref.dat-t(design[batches[[ref]], ] %*% B.hat))^2) %*% rep(1/n.batches[ref],n.batches[ref]) # FIXME
    } else {
      var.pooled <- ((dat-t(design %*% B.hat))^2) %*% rep(1/n.array,n.array) # FIXME
    }
  } else {
    if(!is.null(ref.batch)) {
      ref.dat <- dat[, batches[[ref]]]
      var.pooled <- rowVars(ref.dat-t(design[batches[[ref]], ]%*%B.hat), na.rm=TRUE)
    } else {
      var.pooled <- rowVars(dat-t(design %*% B.hat), na.rm=TRUE)
    }
  }

  stand.mean <- t(grand.mean) %*% t(rep(1,n.array)) # FIXME
  if(!is.null(design)){
    tmp <- design
    tmp[,c(1:n.batch)] <- 0
    stand.mean <- stand.mean+t(tmp %*% B.hat) #FIXME
  }
  s.data <- (dat-stand.mean)/(sqrt(var.pooled) %*% t(rep(1,n.array))) # FIXME

  ##Get regression batch effect parameters
  message("Fitting L/S model and finding priors")
  batch.design <- design[, 1:n.batch]
  if (!NAs){
    gamma.hat <- solve(crossprod(batch.design), tcrossprod(t(batch.design),
                                                           as.matrix(s.data)))
  } else{
    gamma.hat <- apply(s.data, 1, Beta.NA, batch.design) # FIXME
  }
  delta.hat <- NULL
  for (i in batches){
    if(mean.only==TRUE) {
      delta.hat <- rbind(delta.hat,rep(1,nrow(s.data)))
    } else {
      delta.hat <- rbind(delta.hat, rowVars(s.data[,i], na.rm=TRUE))
    }
  }

  ##Find Priors
  gamma.bar <- rowMeans(gamma.hat)
  t2 <- rowVars(gamma.hat)
  a.prior <- apply(delta.hat, 1, aprior) # FIXME
  b.prior <- apply(delta.hat, 1, bprior) # FIXME

  ## Plot empirical and parametric priors

  if (prior.plots && par.prior) {
    old_pars <- par(no.readonly = TRUE)
    on.exit(par(old_pars))
    par(mfrow=c(2,2))

    ## Top left
    tmp <- density(gamma.hat[1,])
    plot(tmp,  type='l', main=expression(paste("Density Plot of First Batch ",  hat(gamma))))
    xx <- seq(min(tmp$x), max(tmp$x), length=100)
    lines(xx,dnorm(xx,gamma.bar[1],sqrt(t2[1])), col=2)

    ## Top Right
    qqnorm(gamma.hat[1,], main=expression(paste("Normal Q-Q Plot of First Batch ", hat(gamma))))
    qqline(gamma.hat[1,], col=2)

    ## Bottom Left
    tmp <- density(delta.hat[1,])
    xx <- seq(min(tmp$x), max(tmp$x), length=100)
    tmp1 <- list(x=xx, y=dinvgamma(xx, a.prior[1], b.prior[1]))
    plot(tmp, typ="l", ylim=c(0, max(tmp$y, tmp1$y)),
         main=expression(paste("Density Plot of First Batch ", hat(delta))))
    lines(tmp1, col=2)

    ## Bottom Right
    invgam <- 1/qgamma(1-ppoints(ncol(delta.hat)), a.prior[1], b.prior[1])
    qqplot(invgam, delta.hat[1,],
           main=expression(paste("Inverse Gamma Q-Q Plot of First Batch ", hat(delta))),
           ylab="Sample Quantiles", xlab="Theoretical Quantiles")
    lines(c(0, max(invgam)), c(0, max(invgam)), col=2)
  }

  ## Find EB batch adjustments

  gamma.star <- delta.star <- matrix(NA, nrow=n.batch, ncol=nrow(s.data))
  if (par.prior) {
    message("Finding parametric adjustments")
    results <- bplapply(1:n.batch, function(i) {
      if (mean.only) {
        gamma.star <- postmean(gamma.hat[i,], gamma.bar[i], 1, 1, t2[i])
        delta.star <- rep(1, nrow(s.data))
      }
      else {
        temp <- it.sol(s.data[, batches[[i]]], gamma.hat[i, ],
                       delta.hat[i, ], gamma.bar[i], t2[i], a.prior[i],
                       b.prior[i])
        gamma.star <- temp[1, ]
        delta.star <- temp[2, ]
      }
      list(gamma.star=gamma.star, delta.star=delta.star)
    }, BPPARAM = BPPARAM)
    for (i in 1:n.batch) {
      gamma.star[i,] <- results[[i]]$gamma.star
      delta.star[i,] <- results[[i]]$delta.star
    }
  }
  else {
    message("Finding nonparametric adjustments")
    results <- bplapply(1:n.batch, function(i) {
      if (mean.only) {
        delta.hat[i, ] = 1
      }
      temp <- int.eprior(as.matrix(s.data[, batches[[i]]]),
                         gamma.hat[i, ], delta.hat[i, ])
      list(gamma.star=temp[1,], delta.star=temp[2,])
    }, BPPARAM = BPPARAM)
    for (i in 1:n.batch) {
      gamma.star[i,] <- results[[i]]$gamma.star
      delta.star[i,] <- results[[i]]$delta.star
    }
  }

  if(!is.null(ref.batch)){
    gamma.star[ref,] <- 0  ## set reference batch mean equal to 0
    delta.star[ref,] <- 1  ## set reference batch variance equal to 1
  }

  ## Normalize the Data ###
  message("Adjusting the Data\n")

  bayesdata <- s.data
  j <- 1
  for (i in batches){
    bayesdata[,i] <- (bayesdata[,i]-t(batch.design[i,]%*%gamma.star))/(sqrt(delta.star[j,])%*%t(rep(1,n.batches[j]))) # FIXME
    j <- j+1
  }

  bayesdata <- (bayesdata*(sqrt(var.pooled)%*%t(rep(1,n.array))))+stand.mean # FIXME

  ## Do not change ref batch at all in reference version
  if(!is.null(ref.batch)){
    bayesdata[, batches[[ref]]] <- dat[, batches[[ref]]]
  }

  ## MODIFIED: This part was removed since the zero variance genes were never removed
  # ## put genes with 0 variance in any batch back in data
  # if (length(zero.rows) > 0) {
  #   dat.orig[keep.rows, ] <- bayesdata
  #   bayesdata <- dat.orig
  # }

  return(bayesdata)
}

#### helper_seq ####
####  Expand a vector into matrix (columns as the original vector)
vec2mat <- function(vec, n_times){
  return(matrix(rep(vec, n_times), ncol=n_times, byrow=FALSE))
}


####  Monte Carlo integration functions
monte_carlo_int_NB <- function(dat, mu, gamma, phi, gene.subset.n){
  weights <- pos_res <- list()
  for(i in 1:nrow(dat)){
    m <- mu[-i,!is.na(dat[i,])]
    x <- dat[i,!is.na(dat[i,])]
    gamma_sub <- gamma[-i]
    phi_sub <- phi[-i]

    # take a subset of genes to do integration - save time
    if(!is.null(gene.subset.n) & is.numeric(gene.subset.n) & length(gene.subset.n)==1){
      if(i==1){cat(sprintf("Using %s random genes for Monte Carlo integration\n", gene.subset.n))}
      mcint_ind <- sample(1:(nrow(dat)-1), gene.subset.n, replace=FALSE)
      m <- m[mcint_ind, ]; gamma_sub <- gamma_sub[mcint_ind]; phi_sub <- phi_sub[mcint_ind]
      G_sub <- gene.subset.n
    }else{
      if(i==1){cat("Using all genes for Monte Carlo integration; the function runs very slow for large number of genes\n")}
      G_sub <- nrow(dat)-1
    }

    #LH <- sapply(1:G_sub, function(j){sum(log2(dnbinom(x, mu=m[j,], size=1/phi_sub[j])+1))})
    LH <- sapply(1:G_sub, function(j){prod(dnbinom(x, mu=m[j,], size=1/phi_sub[j]))})
    LH[is.nan(LH)]=0;
    if(sum(LH)==0 | is.na(sum(LH))){
      pos_res[[i]] <- c(gamma.star=as.numeric(gamma[i]), phi.star=as.numeric(phi[i]))
    }else{
      pos_res[[i]] <- c(gamma.star=sum(gamma_sub*LH)/sum(LH), phi.star=sum(phi_sub*LH)/sum(LH))
    }

    weights[[i]] <- as.matrix(LH/sum(LH))
  }
  pos_res <- do.call(rbind, pos_res)
  weights <- do.call(cbind, weights)
  res <- list(gamma_star=pos_res[, "gamma.star"], phi_star=pos_res[, "phi.star"], weights=weights)
  return(res)
}


####  Match quantiles
match_quantiles <- function(counts_sub, old_mu, old_phi, new_mu, new_phi){
  new_counts_sub <- matrix(NA, nrow=nrow(counts_sub), ncol=ncol(counts_sub))
  for(a in 1:nrow(counts_sub)){
    for(b in 1:ncol(counts_sub)){
      if(counts_sub[a, b] <= 1){
        new_counts_sub[a,b] <- counts_sub[a, b]
      }else{
        tmp_p <- pnbinom(counts_sub[a, b]-1, mu=old_mu[a, b], size=1/old_phi[a])
        if(abs(tmp_p-1)<1e-4){
          new_counts_sub[a,b] <- counts_sub[a, b]
          # for outlier count, if p==1, will return Inf values -> use original count instead
        }else{
          new_counts_sub[a,b] <- 1+qnbinom(tmp_p, mu=new_mu[a, b], size=1/new_phi[a])
        }
      }
    }
  }
  return(new_counts_sub)
}

#### ComBat_seq ####
#' @importFrom edgeR DGEList estimateGLMCommonDisp estimateGLMTagwiseDisp glmFit glmFit.default getOffset
#' @importFrom stats dnbinom lm pnbinom qnbinom
#' @importFrom utils capture.output
.myComBat_seq <- function(counts, batch, group=NULL, covar_mod=NULL, full_mod=TRUE,
                       shrink=FALSE, shrink.disp=FALSE, gene.subset.n=NULL){
  ########  Preparation  ########
  ## Does not support 1 sample per batch yet
  batch <- as.factor(batch)
  if(any(table(batch)<=1)){
    stop("ComBat-seq doesn't support 1 sample per batch yet")
  }

  ## MODIFIED: Batch effect must be removed in all genes
  # # Remove genes with only 0 counts in any batch
  # keep_lst <- lapply(levels(batch), function(b){
  #   which(apply(counts[, batch==b], 1, function(x){!all(x==0)}))
  # })
  # keep <- Reduce(intersect, keep_lst)
  # rm <- setdiff(1:nrow(counts), keep)
  # countsOri <- counts
  # counts <- counts[keep, ]

  # require bioconductor 3.7, edgeR 3.22.1
  dge_obj <- DGEList(counts=counts)

  ## Prepare characteristics on batches
  n_batch <- nlevels(batch)  # number of batches
  batches_ind <- lapply(1:n_batch, function(i){which(batch==levels(batch)[i])}) # list of samples in each batch
  n_batches <- sapply(batches_ind, length)
  #if(any(n_batches==1)){mean_only=TRUE; cat("Note: one batch has only one sample, setting mean.only=TRUE\n")}
  n_sample <- sum(n_batches)
  cat("Found",n_batch,'batches\n')

  ## Make design matrix
  # batch
  batchmod <- model.matrix(~-1+batch)  # colnames: levels(batch)
  # covariate
  group <- as.factor(group)
  if(full_mod & nlevels(group)>1){
    cat("Using full model in ComBat-seq.\n")
    mod <- model.matrix(~group)
  }else{
    cat("Using null model in ComBat-seq.\n")
    mod <- model.matrix(~1, data=as.data.frame(t(counts)))
  }
  # drop intercept in covariate model
  if(!is.null(covar_mod)){
    if(is.data.frame(covar_mod)){
      covar_mod <- do.call(cbind, lapply(1:ncol(covar_mod), function(i){model.matrix(~covar_mod[,i])}))
    }
    covar_mod <- covar_mod[, !apply(covar_mod, 2, function(x){all(x==1)})]
  }
  # bind with biological condition of interest
  mod <- cbind(mod, covar_mod)
  # combine
  design <- cbind(batchmod, mod)

  ## Check for intercept in covariates, and drop if present
  check <- apply(design, 2, function(x) all(x == 1))
  #if(!is.null(ref)){check[ref]=FALSE} ## except don't throw away the reference batch indicator
  design <- as.matrix(design[,!check])
  cat("Adjusting for",ncol(design)-ncol(batchmod),'covariate(s) or covariate level(s)\n')

  ## Check if the design is confounded
  if(qr(design)$rank<ncol(design)){
    #if(ncol(design)<=(n_batch)){stop("Batch variables are redundant! Remove one or more of the batch variables so they are no longer confounded")}
    if(ncol(design)==(n_batch+1)){stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat-Seq")}
    if(ncol(design)>(n_batch+1)){
      if((qr(design[,-c(1:n_batch)])$rank<ncol(design[,-c(1:n_batch)]))){stop('The covariates are confounded! Please remove one or more of the covariates so the design is not confounded')
      }else{stop("At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat-Seq")}}
  }

  ## Check for missing values in count matrix
  NAs = any(is.na(counts))
  if(NAs){cat(c('Found',sum(is.na(counts)),'Missing Data Values\n'),sep=' ')}


  ########  Estimate gene-wise dispersions within each batch  ########
  cat("Estimating dispersions\n")
  ## Estimate common dispersion within each batch as an initial value
  disp_common <- sapply(1:n_batch, function(i){
    if((n_batches[i] <= ncol(design)-ncol(batchmod)+1) | qr(mod[batches_ind[[i]], ])$rank < ncol(mod)){
      # not enough residual degree of freedom
      return(estimateGLMCommonDisp(counts[, batches_ind[[i]]], design=NULL, subset=nrow(counts)))
    }else{
      return(estimateGLMCommonDisp(counts[, batches_ind[[i]]], design=mod[batches_ind[[i]], ], subset=nrow(counts)))
    }
  })

  ## Estimate gene-wise dispersion within each batch
  genewise_disp_lst <- lapply(1:n_batch, function(j){
    if((n_batches[j] <= ncol(design)-ncol(batchmod)+1) | qr(mod[batches_ind[[j]], ])$rank < ncol(mod)){
      # not enough residual degrees of freedom - use the common dispersion
      return(rep(disp_common[j], nrow(counts)))
    }else{
      return(estimateGLMTagwiseDisp(counts[, batches_ind[[j]]], design=mod[batches_ind[[j]], ],
                                    dispersion=disp_common[j], prior.df=0))
    }
  })
  names(genewise_disp_lst) <- paste0('batch', levels(batch))

  ## construct dispersion matrix
  phi_matrix <- matrix(NA, nrow=nrow(counts), ncol=ncol(counts))
  for(k in 1:n_batch){
    phi_matrix[, batches_ind[[k]]] <- vec2mat(genewise_disp_lst[[k]], n_batches[k])
  }


  ########  Estimate parameters from NB GLM  ########
  cat("Fitting the GLM model\n")
  glm_f <- glmFit(dge_obj, design=design, dispersion=phi_matrix, prior.count=1e-4) #no intercept - nonEstimable; compute offset (library sizes) within function
  alpha_g <- glm_f$coefficients[, 1:n_batch] %*% as.matrix(n_batches/n_sample) #compute intercept as batch-size-weighted average from batches
  new_offset <- t(vec2mat(getOffset(dge_obj), nrow(counts))) +   # original offset - sample (library) size
    vec2mat(alpha_g, ncol(counts))  # new offset - gene background expression # getOffset(dge_obj) is the same as log(dge_obj$samples$lib.size)
  glm_f2 <- glmFit.default(dge_obj$counts, design=design, dispersion=phi_matrix, offset=new_offset, prior.count=1e-4)

  gamma_hat <- glm_f2$coefficients[, 1:n_batch]
  mu_hat <- glm_f2$fitted.values
  phi_hat <- do.call(cbind, genewise_disp_lst)


  ########  In each batch, compute posterior estimation through Monte-Carlo integration  ########
  if(shrink){
    cat("Apply shrinkage - computing posterior estimates for parameters\n")
    mcint_fun <- monte_carlo_int_NB
    monte_carlo_res <- lapply(1:n_batch, function(ii){
      if(ii==1){
        mcres <- mcint_fun(dat=counts[, batches_ind[[ii]]], mu=mu_hat[, batches_ind[[ii]]],
                           gamma=gamma_hat[, ii], phi=phi_hat[, ii], gene.subset.n=gene.subset.n)
      }else{
        invisible(capture.output(mcres <- mcint_fun(dat=counts[, batches_ind[[ii]]], mu=mu_hat[, batches_ind[[ii]]],
                                                    gamma=gamma_hat[, ii], phi=phi_hat[, ii], gene.subset.n=gene.subset.n)))
      }
      return(mcres)
    })
    names(monte_carlo_res) <- paste0('batch', levels(batch))

    gamma_star_mat <- lapply(monte_carlo_res, function(res){res$gamma_star})
    gamma_star_mat <- do.call(cbind, gamma_star_mat)
    phi_star_mat <- lapply(monte_carlo_res, function(res){res$phi_star})
    phi_star_mat <- do.call(cbind, phi_star_mat)

    if(!shrink.disp){
      cat("Apply shrinkage to mean only\n")
      phi_star_mat <- phi_hat
    }
  }else{
    cat("Shrinkage off - using GLM estimates for parameters\n")
    gamma_star_mat <- gamma_hat
    phi_star_mat <- phi_hat
  }


  ########  Obtain adjusted batch-free distribution  ########
  mu_star <- matrix(NA, nrow=nrow(counts), ncol=ncol(counts))
  for(jj in 1:n_batch){
    mu_star[, batches_ind[[jj]]] <- exp(log(mu_hat[, batches_ind[[jj]]])-vec2mat(gamma_star_mat[, jj], n_batches[jj]))
  }
  phi_star <- rowMeans(phi_star_mat)


  ########  Adjust the data  ########
  cat("Adjusting the data\n")
  adjust_counts <- matrix(NA, nrow=nrow(counts), ncol=ncol(counts))
  for(kk in 1:n_batch){
    counts_sub <- counts[, batches_ind[[kk]]]
    old_mu <- mu_hat[, batches_ind[[kk]]]
    old_phi <- phi_hat[, kk]
    new_mu <- mu_star[, batches_ind[[kk]]]
    new_phi <- phi_star
    adjust_counts[, batches_ind[[kk]]] <- match_quantiles(counts_sub=counts_sub,
                                                          old_mu=old_mu, old_phi=old_phi,
                                                          new_mu=new_mu, new_phi=new_phi)
  }

  ## MODIFIED: The adjusted data is returned
  dimnames(adjust_counts) <- dimnames(counts)
  return(adjust_counts)

  ## MODIFIED: The zero count genes were never removed
  # ## Add back genes with only 0 counts in any batch (so that dimensions won't change)
  # adjust_counts_whole <- matrix(NA, nrow=nrow(countsOri), ncol=ncol(countsOri))
  # dimnames(adjust_counts_whole) <- dimnames(countsOri)
  # adjust_counts_whole[keep, ] <- adjust_counts
  # adjust_counts_whole[rm, ] <- countsOri[rm, ]
  # return(adjust_counts_whole)
}
