getDesign=unmarked:::getDesign
name_to_ind=unmarked:::name_to_ind

setMethod("getDesign", "unmarkedFrameOccuMulti",
    function(umf, detformulas, stateformulas, maxOrder, na.rm=TRUE, warn=FALSE,
             newdata=NULL, type="state")
{

  #Format formulas
  #Workaround for parameters fixed at 0
  fixed0 <- stateformulas %in% c("~0","0")
  stateformulas[fixed0] <- "~1"

  stateformulas <- lapply(stateformulas,as.formula)
  detformulas <- lapply(detformulas,as.formula)

  #Generate some indices
  S <- length(umf@ylist) # of species
  if(missing(maxOrder)){
    maxOrder <- S
  }
  z <- expand.grid(rep(list(1:0),S))[,S:1] # z matrix
  colnames(z) <- names(umf@ylist)
  M <- nrow(z) # of possible z states

  # f design matrix
  if(maxOrder == 1){
    dmF <- as.matrix(z)
  } else {
    dmF <- model.matrix(as.formula(paste0("~.^",maxOrder,"-1")),z)
  }
  nF <- ncol(dmF) # of f parameters

  J <- ncol(umf@ylist[[1]]) # max # of samples at a site
  N <- nrow(umf@ylist[[1]]) # of sites

  #Check formulas
  if(length(stateformulas) != nF)
    stop(paste(nF,"formulas are required in stateformulas list"))

  if(length(detformulas) != S)
    stop(paste(S,"formulas are required in detformulas list"))

  if(is.null(siteCovs(umf))) {
    site_covs <- data.frame(placeHolderSite = rep(1, N))
  } else {
    site_covs <- siteCovs(umf)
  }

  if(is.null(obsCovs(umf))) {
    obs_covs <- data.frame(placeHolderObs = rep(1, J*N))
  } else {
    obs_covs <- obsCovs(umf)
  }

  #Add site covs to obs covs if we aren't predicting with newdata
  # Record future column names for obsCovs
  col_names <- c(colnames(obs_covs), colnames(site_covs))

  # add site covariates at observation-level
  obs_covs <- cbind(obs_covs, site_covs[rep(1:N, each = J),])
  colnames(obs_covs) <- col_names

  #Re-format ylist
  index <- 1
  ylong <- lapply(umf@ylist, function(x) {
                   colnames(x) <- 1:J
                   x <- cbind(x,site=1:N,species=index)
                   index <<- index+1
                   x
          })
  ylong <- as.data.frame(do.call(rbind,ylong))
  ylong <- reshape(ylong, idvar=c("site", "species"), varying=list(1:J),
                   v.names="value", direction="long")
  ylong <- reshape(ylong, idvar=c("site","time"), v.names="value",
                    timevar="species", direction="wide")
  ylong <- ylong[order(ylong$site, ylong$time), ]

  #Remove missing values
  if(na.rm){
    naSiteCovs <- which(apply(site_covs, 1, function(x) any(is.na(x))))
    if(length(naSiteCovs>0)){
      stop(paste("Missing site covariates at sites:",
                 paste(naSiteCovs,collapse=", ")))
    }

    naY <- apply(ylong, 1, function(x) any(is.na(x)))
    naCov <- apply(obs_covs, 1, function(x) any(is.na(x)))
    navec <- naY | naCov

    sites_with_missingY <- unique(ylong$site[naY])
    sites_with_missingCov <- unique(ylong$site[naCov])

    ylong <- ylong[!navec,,drop=FALSE]
    obs_covs <- obs_covs[!navec,,drop=FALSE]

    no_data_sites <- which(! 1:N %in% ylong$site)
    if(length(no_data_sites>0)){
      stop(paste("All detections and/or detection covariates are missing at sites:",
                  paste(no_data_sites,collapse=", ")))
    }

    if(sum(naY)>0&warn){
      warning(paste("Missing detections at sites:",
                    paste(sites_with_missingY,collapse=", ")))
    }
    if(sum(naCov)>0&warn){
      warning(paste("Missing detection covariate values at sites:",
                    paste(sites_with_missingCov,collapse=", ")))
    }

  }

  #Start-stop indices for sites
  yStart <- c(1,1+which(diff(ylong$site)!=0))
  yStop <- c(yStart[2:length(yStart)]-1,nrow(ylong))

  y <- as.matrix(ylong[,3:ncol(ylong)])

  #Indicator matrix for no detections at a site
  Iy0 <- do.call(cbind, lapply(umf@ylist,
                               function(x) as.numeric(rowSums(x, na.rm=T)==0)))

  #Save formatted covariate frames for use in model frames
  #For predicting with formulas etc
  site_ref <- site_covs
  obs_ref <- obs_covs

  #Assign newdata as the covariate frame if it is provided
  if(!is.null(newdata)){
    if(type == "state"){
      site_covs <- newdata
    } else if(type == "det"){
      obs_covs <- newdata
    }
  }

  #Design matrices + parameter counts
  #For f/occupancy
  fInd <- c()
  sf_no0 <- stateformulas[!fixed0]
  var_names <- colnames(dmF)[!fixed0]
  dmOcc <- lapply(seq_along(sf_no0),function(i){
                    fac_col <- site_ref[, sapply(site_ref, is.factor), drop=FALSE]
                    mf <- model.frame(sf_no0[[i]], site_ref)
                    xlevs <- lapply(fac_col, levels)
                    xlevs <- xlevs[names(xlevs) %in% names(mf)]
                    out <- model.matrix(sf_no0[[i]],
                                        model.frame(stats::terms(mf), site_covs, na.action=stats::na.pass, xlev=xlevs))
                    colnames(out) <- paste('[',var_names[i],'] ',
                                           colnames(out), sep='')
                    fInd <<- c(fInd,rep(i,ncol(out)))
                    out
          })
  fStart <- c(1,1+which(diff(fInd)!=0))
  fStop <- c(fStart[2:length(fStart)]-1,length(fInd))
  occParams <- unlist(lapply(dmOcc,colnames))
  nOP <- length(occParams)

  #For detection
  dInd <- c()
  dmDet <- lapply(seq_along(detformulas),function(i){
                    fac_col <- obs_ref[, sapply(obs_ref, is.factor), drop=FALSE]
                    mf <- model.frame(detformulas[[i]], obs_ref)
                    xlevs <- lapply(fac_col, levels)
                    xlevs <- xlevs[names(xlevs) %in% names(mf)]
                    out <- model.matrix(detformulas[[i]],
                                        model.frame(stats::terms(mf), obs_covs, na.action=stats::na.pass, xlev=xlevs))
                    colnames(out) <- paste('[',names(umf@ylist)[i],'] ',
                                           colnames(out),sep='')
                    dInd <<- c(dInd,rep(i,ncol(out)))
                    out
          })
  dStart <- c(1,1+which(diff(dInd)!=0)) + nOP
  dStop <- c(dStart[2:length(dStart)]-1,length(dInd)+nOP)
  detParams <- unlist(lapply(dmDet,colnames))
  #nD <- length(detParams)

  #Combined
  paramNames <- c(occParams,detParams)
  nP <- length(paramNames)

  mget(c("N","S","J","M","nF","fStart","fStop","fixed0","dmF","dmOcc","dmDet",
         "dStart","dStop","y","yStart","yStop","Iy0","z","nOP","nP","paramNames"))
})




setMethod("predict", "unmarkedFitOccuMulti",
     function(object, type, newdata,
              #backTransform = TRUE, na.rm = TRUE,
              #appendData = FALSE,
              se.fit=TRUE, level=0.95, species=NULL, cond=NULL, nsims=100,
              ...)
  {

  type <- match.arg(type, c("state", "det"))

  if(is.null(hessian(object))){
    se.fit = FALSE
  }

  species <- name_to_ind(species, names(object@data@ylist))
  cond <- name_to_ind(cond, names(object@data@ylist))

  if(missing(newdata)){
    newdata <- NULL
  } else {
    if(! class(newdata) %in% c('data.frame')){
      stop("newdata must be a data frame")
    }
  }

  maxOrder <- object@call$maxOrder
  if(is.null(maxOrder)) maxOrder <- length(object@data@ylist)

  dm <- getDesign(object@data,object@detformulas,object@stateformulas,
                  maxOrder, na.rm=F, newdata=newdata, type=type)

  params <- coef(object)
  low_bound <- (1-level)/2
  up_bound <- level + (1-level)/2


  if(type=="state"){
    N <- nrow(dm$dmOcc[[1]]); nF <- dm$nF; dmOcc <- dm$dmOcc;
    fStart <- dm$fStart; fStop <- dm$fStop; fixed0 <- dm$fixed0
    t_dmF <- t(dm$dmF)

    calc_psi <- function(params){

      f <- matrix(NA,nrow=N,ncol=nF)
      index <- 1
      for (i in 1:nF){
        if(fixed0[i]){
          f[,i] <- 0
        } else {
          f[,i] <- dmOcc[[index]] %*% params[fStart[index]:fStop[index]]
          index <- index + 1
        }
      }
      psi <- exp(f %*% t_dmF)
      as.matrix(psi/rowSums(psi))
    }

    psi_est <- calc_psi(params)

    if(se.fit){
      cat('Bootstrapping confidence intervals with',nsims,'samples\n')
      Sigma <- vcov(object)
      samp <- array(NA,c(dim(psi_est),nsims))
      for (i in 1:nsims){
        samp[,,i] <- calc_psi(MASS::mvrnorm(1, coef(object), Sigma))
      }
    }

    if(!is.null(species)){

      sel_col <- species

      if(!is.null(cond)){
        if(any(sel_col %in% abs(cond))){
          stop("Species can't be conditional on itself")
        }
        ftemp <- object@data@fDesign
        swap <- -1*cond[which(cond<0)]
        ftemp[,swap] <- 1 - ftemp[,swap]
        num_inds <- apply(ftemp[,c(sel_col,abs(cond))] == 1,1,all)
        denom_inds <- apply(ftemp[,abs(cond),drop=F] == 1,1,all)
        est <- rowSums(psi_est[,num_inds,drop=F]) /
          rowSums(psi_est[,denom_inds, drop=F])
        if(se.fit){
          samp_num <- apply(samp[,num_inds,,drop=F],3,rowSums)
          samp_denom <- apply(samp[,denom_inds,,drop=F],3,rowSums)
          samp <- samp_num / samp_denom
        }

      } else {
        num_inds <- apply(object@data@fDesign[,sel_col,drop=FALSE] == 1,1,all)
        est <- rowSums(psi_est[,num_inds,drop=F])
        if(se.fit){
          samp <- samp[,num_inds,,drop=F]
          samp <- apply(samp, 3, rowSums)
        }
      }

      if(se.fit){
        if(!is.matrix(samp)) samp <- matrix(samp, nrow=1)
        boot_se <- apply(samp,1,sd, na.rm=T)
        boot_low <- apply(samp,1,quantile,low_bound, na.rm=T)
        boot_up <- apply(samp,1,quantile,up_bound, na.rm=T)
      } else{
        boot_se <- boot_low <- boot_up <- NA
      }
      return(data.frame(Predicted=est,
                        SE=boot_se,
                        lower=boot_low,
                        upper=boot_up))

    } else {
      codes <- apply(dm$z,1,function(x) paste(x,collapse=""))
      colnames(psi_est)  <- paste('psi[',codes,']',sep='')
      if(se.fit){
        boot_se <- apply(samp,c(1,2),sd, na.rm=T)
        boot_low <- apply(samp,c(1,2),quantile,low_bound, na.rm=T)
        boot_up <- apply(samp,c(1,2),quantile,up_bound, na.rm=T)
        colnames(boot_se) <- colnames(boot_low) <- colnames(boot_up) <-
          colnames(psi_est)
      } else {
        boot_se <- boot_low <- boot_up <- NA
      }
      return(list(Predicted=psi_est,
                  SE=boot_se,
                  lower=boot_low,
                  upper=boot_up))
    }
  }

  if(type=="det"){
    #based on
    #https://blog.methodsconsultants.com/posts/delta-method-standard-errors/
    S <- dm$S; dmDet <- dm$dmDet
    dStart <- dm$dStart; dStop <- dm$dStop

    out <- list()
    z <- qnorm(low_bound,lower.tail=F)
    N <- nrow(dmDet[[1]])
    for (i in 1:S){

      inds <- dStart[i]:dStop[i]
      param_sub <- params[inds]
      est <- plogis(dmDet[[i]] %*% param_sub)

      if(se.fit){
        cov_sub <- vcov(object)[inds-sum(dm$fixed0),inds-sum(dm$fixed0)]
        se_est <- lower <- upper <- numeric(N)
        for (j in 1:N){
          x <- dmDet[[i]][j,]
          xb <- stats::dlogis(t(x) %*% param_sub)
          v <- xb %*% t(x) %*% cov_sub %*% x %*% xb
          se_est[j] <- sqrt(v)
          lower[j] <- est[j] - z*se_est[j]
          upper[j] <- est[j] + z*se_est[j]
        }
      } else {
        se_est <- lower <- upper <- NA
      }
      out[[i]] <- data.frame(Predicted=est,SE=se_est,
                            lower=lower,upper=upper)
    }
    names(out) <- names(object@data@ylist)
    if(!is.null(species)){
      return(out[[species]])
    }
    return(out)
  }
  stop("type must be 'det' or 'state'")
})
