
require(clinfun);

##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
##                              TOOLS
##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
sub.txt.w.no <- function(numbers, template.f, out.f="rst.txt", sub.str="AA") {
    if (!file.exists(template.f)) {
        return;
    }
    ##read
    tpla <- readChar(template.f, file.info(template.f)$size);

    ##subsitute
    for (i in 1:length(numbers)) {
        tpla <- sub(sub.str, numbers[i], tpla);
    }

    ##write out
    write(tpla, file=out.f);
}


##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
##                       BATCH EFFECT
##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
get.truth <- function(mu, par.error, p0, p1, beta = 0, n.total = 100000) {
    true.pts <- simu.all.pts(batch.size = 1,
                             n.total = n.total,
                             mu = mu,
                             beta = beta,
                             par.error = par.error);
    ## measure 1: ratio
    var.all     <- var(true.pts$Y);
    ## variance without batch effect
    var.nobatch <- var(true.pts$Y.NoBatch);
    var.ratio   <- var.all/var.nobatch - 1;

    ## cut off
    cut.y <- NULL;
    for (p in c(p0,p1)) {
        cut.y <- c(cut.y, as.numeric(quantile(true.pts$Y, 1-p)));
    }

    c(var.ratio   = var.ratio,
      cuty.h0     = cut.y[1],
      cuty.h1     = cut.y[2],
      var.nobatch = var.nobatch,
      var.all     = var.all);
}

##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
##                       SIMULATION
##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------

## simulate random effects gamma, delta and epsilon
simu.error <- function(n.simu,
                       error.type = c("normal", "skewed"),
                       ysig   = 1,
                       mu     = 0,
                       skew.n = NULL, skew.p = NULL, skew.noise = 0.0001) {

    type <- match.arg(error.type);
    rst <- switch(type,
                  normal = {rnorm(n.simu, mu, ysig)},
                  skewed = {
        smu       <- skew.n * (1-skew.p) / skew.p;
        va        <- skew.n * (1-skew.p) / skew.p^2;
        noise.sig <- skew.noise;
        rst       <- rnbinom(n = n.simu, size = skew.n, prob = skew.p);
        rst       <- rst - smu + rnorm(n.simu, mu, noise.sig);
        rst       <- rst/sqrt(va + noise.sig^2)*ysig;
    });

    rst
}

simu.batch.error <- function(batch.size, n.batch,
                             par.error  = list(gamma   = list(error.type = "normal", ysig = 1),
                                               delta   = list(error.type = "normal", ysig = 1, mu = 1),
                                               epsilon = list(error.type = "normal", ysig = 1))) {

    beffs <- NULL;
    for (ef in c("gamma", "delta")) {
        cur.par <- par.error[[ef]];
        cur.eff <- do.call(simu.error, c(n.simu = n.batch, cur.par));
        beffs   <- cbind(beffs, cur.eff);
    }

    gamma   <- beffs[,1];
    delta   <- beffs[,2];
    epsilon <- do.call(simu.error,
                       c(n.simu = batch.size * n.batch, par.error[["epsilon"]]));

    error <- rep(gamma, each = batch.size) + rep(delta, each = batch.size) * epsilon;

    list(gamma   = gamma,
         delta   = delta,
         epsilon = epsilon,
         error   = error);
}


## simulate for an individual batch
simu.all.pts <- function(n.total, batch.size = 1, n.batch = NULL, mu = c(0, 0), beta = 0, ...) {

    stopifnot(n.total > 0);

    if (!xor(is.null(batch.size), is.null(n.batch)))
        stop("Please specify either batch.size or n.batch.")

    if (is.null(batch.size)) {
        batch.size <- ceiling(n.total/n.batch);
    } else {
        n.batch <- as.numeric(ceiling(n.total / batch.size));
    }

    nt             <- length(mu);
    ystar          <- NULL;
    ystar.nobatch  <- NULL;
    gamma          <- NULL;
    delta          <- NULL;
    epsilon        <- NULL;
    for (i in 1:nt) {
        cur.mu <- mu[i];
        if (i > 1)
            cur.mu <- cur.mu + ystar.nobatch[, i-1];

        cur.error <- simu.batch.error(batch.size, n.batch, ...);

        gamma    <- cbind(gamma,   cur.error$gamma);
        delta    <- cbind(delta,   cur.error$delta);
        epsiolon <- cbind(epsilon, cur.error$epsilon);

        ystar         <- cbind(ystar,         cur.mu + cur.error$error);
        ystar.nobatch <- cbind(ystar.nobatch, cur.mu + cur.error$epsilon);
    }

    if (1 == nt) {
        ##absolute outcome
        Y         <- exp(ystar);
        Y.NoBatch <- exp(ystar.nobatch);
    } else {
        ## ratio
        Y <- exp(ystar[,nt])/exp(ystar[,1]);
        Y.NoBatch <- exp(ystar.nobatch[,nt])/exp(ystar.nobatch[,1]);
    }

    list(Y          = Y[1:n.total],
         Y.NoBatch  = Y.NoBatch[1:n.total],
         Ystar      = ystar[1:n.total,],
         Ystar.nob  = ystar.nobatch[1:n.total,],
         mu         = mu,
         beta       = beta,
         batch      = rep(1:n.batch, each = batch.size)[1:n.total],
         batch.size = batch.size,
         n.batch    = n.batch,
         gamma      = gamma,
         delta      = delta,
         epsilon    = epsilon);
}


## simulate trials
simu.trial <- function(p1 = 0.3, p0 = 0.15,
                       alpha = 0.05, power = 0.8,
                       batch.size = 3, n.batch = NULL,
                       nreps = 10000,
                       mu = 0, beta = 0,
                       par.error = list(gamma   = list(error.type = "normal", ysig = 1),
                                        delta   = list(error.type = "normal", ysig = 1, mu = 1),
                                        epsilon = list(error.type = "normal", ysig = 1)),
                       cuty.h01 = NULL, n.truepts = 1000000, seed = NULL) {

    if (is.numeric(seed))
        set.seed(seed);

    ## truth
    truth <- get.truth(mu, par.error, p0, p1, beta = beta, n.total = n.truepts);

    if (is.null(cuty.h01))
        cuty.h01 <- truth[c("cuty.h0", "cuty.h1")];

    ## designs
    design.simon  <- get.simon.2.opt(pu = p0, pa = p1, ep1   = alpha, ep2  = 1 - power);
    design.single <- get.single.size(p0 = p0, p1 = p1, alpha = alpha, beta = 1 - power);

    ## replications
    all.rst <- NULL;
    for (rep in 1:NREPS) {
        cur.rst <- NULL;

        ##---- simon design ------------
        simon.n      <- c(design.simon["n1"],
                          design.simon["n"] - design.simon["n1"]);
        cur.simon.pt <- NULL;
        for (sn in simon.n) {
            cur.pt <- simu.all.pts(n.total    = sn,
                                   batch.size = batch.size,
                                   n.batch    = n.batch,
                                   mu         = mu,
                                   beta       = beta,
                                   par.error  = par.error);
            cur.simon.pt <- c(cur.simon.pt, cur.pt$Y);
        }

        for (cuty in cuty.h01) {
            cur.rst <- c(cur.rst,
                         sum.simon(cur.simon.pt > cuty, design.simon, alpha, p0));
        }

        ##----single stage design -----
        cur.single.pt <- simu.all.pts(n.total    = as.numeric(design.single["n"]),
                                      batch.size = batch.size,
                                      n.batch    = n.batch,
                                      mu         = mu,
                                      beta       = beta,
                                      par.error  = par.error);


        for (cuty in cuty.h01) {
            cur.rst <- c(cur.rst,
                         sum.single(cur.single.pt$Y > cuty, alpha, p0));
        }

        all.rst <- rbind(all.rst, cur.rst);
    }

    rst <- list(var.ratio     = truth["var.ratio"],
                sum.simon.h0  = sum.simon.all(all.rst[,1:6],   truep = p0),
                sum.simon.h1  = sum.simon.all(all.rst[,7:12],  truep = p1),
                sum.single.h0 = sum.single.all(all.rst[,13:17],truep = p0),
                sum.single.h1 = sum.single.all(all.rst[,18:22],truep = p1),
                design.simon  = design.simon,
                design.single = design.single,
                p01           = c(p0, p1),
                mu            = mu,
                batch.size    = batch.size,
                n.batch       = n.batch,
                par.error     = par.error)
}

##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
##                       SINGLE-STAGE
##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------

get.z.ab <- function(alpha=0.05, beta=0.2) {
    z.alpha <- qnorm(1 - alpha/2);
    z.beta  <- qnorm(1 - beta);

    c(z.alpha, z.beta);
}


get.single.size <- function(p1, p0, alpha=0.05, beta=0.2) {
    zab     <- get.z.ab(alpha=alpha, beta=beta);
    z.alpha <- zab[1];
    z.beta  <- zab[2];
    N       <- (z.alpha + z.beta)^2/(p1 - p0)^2;
    N       <- N * p1 * (1-p1);

    ##return
    c(n = ceiling(N));
}

sum.single <- function(resp, alpha, p0) {
    cur.test <- binom.test(sum(resp), length(resp), p = p0, conf.level = 1-alpha);
    rej      <- cur.test$conf.int[1] > p0 | cur.test$conf.int[2] < p0;

    c(cur.test$estimate,
      cur.test$p.value < alpha,
      cur.test$conf.int,
      rej);
}

sum.single.all <- function(single.rst, truep) {
    all.mean   <- as.numeric(apply(single.rst, 2, mean));
    bias       <- all.mean[1] - truep;
    mse        <- mean((single.rst[,1] - truep)^2);
    conf.width <- mean(single.rst[,4] - single.rst[,3]);

    c(AvgCIW  = conf.width,
      RejRate = all.mean[5],
      Bias    = bias,
      MSE     = mse,
      RejPval = all.mean[2]);
}


##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
##                       SIMON 2-STAGE
##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------

##simon two stage optimal design
get.simon.2.opt <- function(...) {
    ss    <- ph2simon(...);
    nopt  <- which.min(ss$out[,5]);
    ss$out[nopt,];
}

##simons two stage
get.simon.2stage <- function(ntotal, p0, p1, alpha=0.025, beta=0.2) {
    f.pr <- function(n1,r1,n,r,p) {
        rst <- pbinom(r1, n1, p);
        for (i in (r1+1):min(n1,r)) {
            rst <- rst + dbinom(i, n1, p) * pbinom(r - i, n-n1, p);
        }
        1-rst
    }

    rst <- NULL;
    for (n in ntotal) {
        for (n1 in 1:(n-1)) {
            for (r in 0:n) {
                cat("n1 =  ", n1, ", r =", r, "\n");
                for (r1 in 0:min(r,n1)) {
                    p.alpha <- f.pr(n1,r1,n,r,p0);
                    p.power <- f.pr(n1,r1,n,r,p1);
                    rst     <- rbind(rst, c(n1,r1,n,r,p.alpha, p.power));
                }
            }
        }
    }

    ## constraints
    inx <- which(rst[,5] < alpha & rst[,6] > 1 - beta);
    rst <- rst[inx,];

    ## early stopping
    pet <- pbinom(rst[,2], rst[,1], p0);
    en  <- pet * rst[,1] + (1-pet) * rst[,3];
    rst <- cbind(rst, pet, en);

    colnames(rst) <- c("N1", "R1", "N", "R", "Alpha", "Power", "EarlyStopping", "EN");
    rst
}

sum.simon <- function(resp, design, alpha, p0) {
    stop.1 <- sum(resp[1:design["n1"]]) <= design["r1"];
    rej    <- sum(resp) > design["r"];
    enroll <- ifelse(stop.1, design["n1"], design["n"]);

    resp.obj <- ifelse(stop.1, resp[1:design["n1"]], resp);
    cur.test <- binom.test(sum(resp.obj), length(resp.obj),
                           p = p0,
                           conf.level = 1-alpha);

    c(stop1  = stop.1,
      rej    = rej,
      enroll = enroll,
      cur.test$estimate,
      cur.test$conf.int);
}

sum.simon.all <- function(simon.rst, truep) {
    all.mean <- as.numeric(apply(simon.rst, 2, mean));
    bias     <- all.mean[4] - truep;
    mse      <- mean((simon.rst[,4] - truep)^2);
    ciw      <- mean(simon.rst[,6] - simon.rst[,5]);

    c(AvgN = all.mean[3],
      AvgEarlyStop = all.mean[1],
      AvgCIW = ciw,
      RejRate = all.mean[2],
      Bias = bias,
      MSE = mse);
}

##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
##                       SIMULATION RESULTS
##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
simu.combine <- function(path    = "./Results/absoutcome/",
                         sce     = 1:4,
                         batch   = c(paste("batchsize", c(1,3,6,9), sep=""), "nbatch2", "nbatch1"),
                         vratio  = c("0.05", "0.10", "0.15", "0.20"),
                         postfix = c("", "_gamma_fixed", "_gamma_skewed")) {
    overall.rst <- NULL;
    for (s in sce) {
        for (bs in batch) {
            for (v in vratio) {
                for (pf in postfix) {
                    ## load data
                    cur.f <- paste(path, "simu_sce", s, "_", bs, "_vratio", v, pf, ".Rdata", sep="");
                    load(cur.f);

                    if (!is.null(rst$batch.size)) {
                        bsize <- rst$batch.size;
                    } else if (2 == rst$n.batch) {
                        bsize <- "B2";
                    } else {
                        bsize <- "BF";
                    }

                    ## organize results
                    cur.rst     <- data.frame(VRatio = v,
                                              Sce    = s,
                                              Study  = pf,
                                              p0     = rst$p01[1],
                                              p1     = rst$p01[2],
                                              Batch  = bsize,
                                              MSize  = as.numeric(rst$design.simon["n"]),
                                              MExp   = as.numeric(rst$design.simon["EN(p0)"]),
                                              M0N    = as.numeric(rst$sum.simon.h0["AvgN"]),
                                              M0Rej  = as.numeric(rst$sum.simon.h0["RejRate"]),
                                              M1N    = as.numeric(rst$sum.simon.h1["AvgN"]),
                                              M1Rej  = as.numeric(rst$sum.simon.h1["RejRate"]),
                                              S0CIW  = as.numeric(rst$sum.single.h0["AvgCIW"]),
                                              S0Rej  = as.numeric(rst$sum.single.h0["RejRate"]),
                                              S0MSE  = as.numeric(rst$sum.single.h0["MSE"]),
                                              S0Bias = as.numeric(rst$sum.single.h0["Bias"]),
                                              S1CIW  = as.numeric(rst$sum.single.h1["AvgCIW"]),
                                              S1Rej  = as.numeric(rst$sum.single.h1["RejRate"]),
                                              S1MSE  = as.numeric(rst$sum.single.h1["MSE"]),
                                              S1Bias = as.numeric(rst$sum.single.h1["Bias"]));

                    ## append
                    overall.rst <- rbind(overall.rst, cur.rst);
                }
            }
        }
    }

    row.names(overall.rst) <- NULL;
    data.frame(overall.rst);
}

plot.simu <- function(rst.data, y, ylab, yinter = NULL,  ylim = NULL) {

    rst <- ggplot(data = rst.data, aes(x = Batch, group = VRatio)) +
        geom_line(aes_string(y = y, col = "VRatio")) +
        theme_bw() +
        labs(y = ylab, x = "Batch Size") +
        theme(legend.position = "top") +
        facet_grid(Group~Sce);

    if (!is.null(ylim))
        rst <- rst + scale_y_continuous(limits = ylim);

    if (!is.null(yinter))
        rst <- rst + geom_hline(yintercept = yinter, linetype="dashed", col = "red");


    rst
}
