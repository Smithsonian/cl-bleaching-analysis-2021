#################################################################
##  Modification of function ME() from the spatialreg package  ##
#################################################################
# There are two changes here from the ME() function
# in spatialreg [Bivand, R., and Piras, G. (2015). Comparing implementations of estimation methods for spatial econometrics. Journal of Statistical Software 63: 1-36]
# One is that we include a maximum number of eigenvectors that can be
# added to the glm model, by adding a counter in the while-loop below.
# This is strictly a time- and computation-saving change, since our
# model selection criterion is to find, across all parameterizations of the
# spatial weights matrix under consideration, the most parsimonious model
# (i.e., the model with fewest eigenvectors) that yields a non-significant
# Moran's I statistic, we only need to add enough eigenvectors so that at
# least one parameter combination meets this criterion.

# The second change is that we use the MEM function, rather than the eigen()
# function, to extract only the positive eigenvectors. We do this following
# Bauman et al 2018 (Ecology 99: 2159-2166, see MS for full citation), since
# on ecological grounds we expect autocorrelation to be positive.

# The last change is that we save Moran's I, rather than the z-score, from
# each fitted GLM.

# You can modify the spatialreg ME() function by using the function call 
# trace(ME, edit = T) and then replacing everything in that window with 
# the code below.
######################################################################## 

function (formula, data = list(), family = gaussian, weights, 
          offset, na.action = na.fail, listw = NULL, alpha = 0.05, 
          nsim = 99, verbose = NULL, stdev = FALSE, zero.policy = NULL) 
{
        MoraneI.boot <- function(var, i, ...) {
                var <- var[i]
                I <- (n/S0) * (crossprod(sW %*% var, var))/cpvar
                return(c(as(I, "matrix")))
        }
        MIR_a <- function(resids, sW, n, cpvar, S0, nsim, stdev = TRUE, 
                          par_boot_args = list()) {
                boot1 <- boot(resids, statistic = MoraneI.boot, R = nsim, 
                              sim = "permutation", sW = sW, n = n, S0 = S0, cpvar = cpvar, 
                              parallel = par_boot_args$parallel, ncpus = par_boot_args$ncpus, 
                              cl = par_boot_args$cl)
                mi <- boot1$t0
                if (stdev) {
                        zi <- (boot1$t0 - mean(boot1$t))/sqrt(var(boot1$t))
                        pri <- pnorm(abs(zi), lower.tail = FALSE)
                }
                else {
                        zi <- NA
                        pri <- (sum(boot1$t >= mi) + 1)/(nsim + 1)
                }
                res <- list(estimate = mi, statistic = zi, p.value = pri)
                res
        }
        if (is.null(zero.policy)) 
                zero.policy <- get("zeroPolicy", envir = .spatialregOptions)
        stopifnot(is.logical(zero.policy))
        if (is.null(verbose)) 
                verbose <- get("verbose", envir = .spatialregOptions)
        stopifnot(is.logical(verbose))
        if (is.character(family)) 
                family <- get(family, mode = "function", envir = parent.frame())
        if (is.function(family)) 
                family <- family()
        if (is.null(family$family)) {
                print(family)
                stop("'family' not recognized")
        }
        mf <- match.call(expand.dots = FALSE)
        m <- match(c("formula", "data", "weights", "offset", "na.action"), 
                   names(mf), 0)
        mf <- mf[c(1, m)]
        mf$drop.unused.levels <- TRUE
        mf[[1]] <- quote(stats::model.frame)
        mf <- eval(mf, parent.frame())
        mt <- attr(mf, "terms")
        Y <- model.extract(mf, "response")
        X <- model.matrix(mt, mf)
        weights <- model.weights(mf)
        if (!is.null(weights) && any(weights < 0)) 
                stop("negative weights not allowed")
        offset <- model.offset(mf)
        if (!is.null(offset) && length(offset) != NROW(Y)) 
                stop("number of offsets should equal number of observations")
        na.act <- attr(mf, "na.action")
        if (is.null(listw)) 
                stop("listw required")
        if (!is.null(na.act)) {
                subset <- !(1:length(listw$neighbours) %in% na.act)
                listw <- subset(listw, subset, zero.policy = zero.policy)
        }
        listw <- listw2U(listw)
        sW <- as(listw, "CsparseMatrix")
        Wmat <- listw2mat(listw)
        n <- ncol(Wmat)
        S0 <- Szero(listw)
        glm_fit <- glm.fit(x = X, y = Y, weights = weights, offset = offset, 
                           family = family)
        glm_res <- glm_fit$y - glm_fit$fitted.values
        cpvar <- crossprod(glm_res)
        cores <- get.coresOption()
        if (is.null(cores)) {
                parallel <- "no"
        }
        else {
                parallel <- ifelse(get.mcOption(), "multicore", "snow")
        }
        ncpus <- ifelse(is.null(cores), 1L, cores)
        cl <- NULL
        if (parallel == "snow") {
                cl <- get.ClusterOption()
                if (is.null(cl)) {
                        parallel <- "no"
                        warning("no cluster in ClusterOption, parallel set to no")
                }
        }
        par_boot_args <- list(parallel = parallel, ncpus = ncpus, 
                              cl = cl)
        mRES <- MIR_a(glm_res, sW = sW, n = n, cpvar = cpvar, S0 = S0, 
                      nsim = nsim, stdev = stdev, par_boot_args = par_boot_args)
        pIZ <- mRES$p.value
        tres <- c(NA, mRES$estimate, pIZ)
        if (pIZ > alpha) 
                stop("base correlation larger than alpha")
        Cent <- diag(n) - (matrix(1/n, n, n))
        CWC <- Cent %*% Wmat %*% Cent
        rm(Cent, Wmat)
        CWC2 <- 0.5 * (CWC + t(CWC))
        rm(CWC)
        # Here instead of using the eigen() function, we use the mem() function of the adespatial package to
        # calculate the "positive" eigenvectors.
        eV <- mem(listw = listw, MEM.autocor = "positive")
        rm(CWC2)
        iZ <- numeric(ncol(eV))
        for (i in 1:ncol(eV)) {
                iX <- cbind(X, eV[, i])
                i_glm <- glm.fit(x = iX, y = Y, weights = weights, offset = offset, 
                                 family = family)
                glm_res <- i_glm$y - i_glm$fitted.values
                cpvar <- crossprod(glm_res)
                iZ[i] <- c(as((n/S0) * (crossprod(sW %*% glm_res, glm_res))/cpvar, 
                              "matrix"))
        }
        min_iZ <- which.min(abs(iZ))
        X <- cbind(X, eV[, min_iZ])
        glm_fit <- glm.fit(x = X, y = Y, weights = weights, offset = offset, 
                           family = family)
        glm_res <- glm_fit$y - glm_fit$fitted.values
        cpvar <- crossprod(glm_res)
        mRES <- MIR_a(glm_res, sW = sW, n = n, cpvar = cpvar, S0 = S0, 
                      nsim = nsim, stdev = stdev, par_boot_args = par_boot_args)
        pIZ <- mRES$p.value
        used <- rep(FALSE, n)
        used[min_iZ] <- TRUE
        min_v <- min_iZ
        if (verbose) 
                cat("eV[,", min_iZ, "], I: ", mRES$estimate, " ZI: ", 
                    mRES$statistic, ", pr(ZI): ", pIZ, "\n", sep = "")
        tres <- rbind(tres, c(min_iZ, mRES$estimate, pIZ))
        # We add a counter that also controls the while-loop. If the number of eigenvectors goes
        # above a certain threshold then it should stop looking for more. 
        counter <- 0
        while ((pIZ <= alpha) & (counter <=8)) { 
                for (i in 1:ncol(eV)) {
                        if (!used[i]) {
                                iX <- cbind(X, eV[, i])
                                i_glm <- glm.fit(x = iX, y = Y, weights = weights, 
                                                 offset = offset, family = family)
                                glm_res <- i_glm$y - i_glm$fitted.values
                                cpvar <- crossprod(glm_res)
                                iZ[i] <- c(as((n/S0) * (crossprod(sW %*% glm_res, 
                                                                  glm_res))/cpvar, "matrix"))
                        }
                        else iZ[i] <- NA
                }
                min_iZ <- which.min(abs(iZ))
                X <- cbind(X, eV[, min_iZ])
                glm_fit <- glm.fit(x = X, y = Y, weights = weights, 
                                   offset = offset, family = family)
                glm_res <- glm_fit$y - glm_fit$fitted.values
                cpvar <- crossprod(glm_res)
                mRES <- MIR_a(glm_res, sW = sW, n = n, cpvar = cpvar, 
                              S0 = S0, nsim = nsim, stdev = stdev, par_boot_args = par_boot_args)
                pIZ <- mRES$p.value
                used[min_iZ] <- TRUE
                min_v <- c(min_v, min_iZ)
                if (verbose) 
                        cat("eV[,", min_iZ, "], I: ", mRES$estimate, " ZI: ", 
                            mRES$statistic, ", pr(ZI): ", pIZ, "\n", sep = "")
                # Instead of storing the Moran's I z-score in the results, we store the
                # Moran's I value. This is referring to mRES$estimate.
                tres <- rbind(tres, c(min_iZ, mRES$estimate, pIZ))
                counter <- counter + 1 # We increase the counter in each iteration.
        }
        sv <- eV[, min_v, drop = FALSE]
        colnames(sv) <- paste("vec", min_v, sep = "")
        colnames(tres) <- c("Eigenvector", "moranI", "pr(ZI)")
        rownames(tres) <- 0:(nrow(tres) - 1)
        res <- list(selection = tres, vectors = sv)
        if (!is.null(na.act)) 
                res$na.action <- na.act
        class(res) <- "Me_res"
        res
}
