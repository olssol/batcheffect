
require(clinfun);

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
        mu        <- skew.n * (1-skew.p) / skew.p;
        va        <- skew.n * (1-skew.p) / skew.p^2;
        noise.sig <- skew.noise;
        rst       <- rnbinom(n.simu, skew.n, skew.p);
        rst       <- rst - mu + rnorm(n.simu, mu, noise.sig);
        rst       <- rst/sqrt(va + noise.sig^2)*ysig;
    });

    rst
}

## simulate for an individual batch
simu.batch <- function(batch.size, mu = 0,
                       par.error = list(gamma   = list(error.type = "normal", ysig = 1),
                                        delta   = list(error.type = "normal", ysig = 1, mu = 1),
                                        epsilon = list(error.type = "normal", ysig = 1)),
                       take.exp = TRUE) {
    beffs <- NULL;
    for (ef in c("gamma", "delta")) {
        cur.par <- par.error[[ef]];
        cur.eff <- do.call(simu.error, c(n.simu = 1, cur.par));
        beffs   <- c(beffs, cur.eff);
    }

    epsilons <- do.call(simu.error,
                        c(n.simu = batch.size, par.error[["epsilon"]]));
    rst <- mu + beffs[1] + beffs[2] * epsilons;

    if (take.exp)
        rst <- exp(rst);

    rst
}

## simu all patients in batch size B
simu.all.pts <- function(n.tot, batch.size, ...) {
    rst   <- NULL;
    cur.n <- 0;
    cur.b <- 0;
    while (cur.n < n.tot) {
        cur.b   <- cur.b + 1;
        cur.rst <- simu.batch(batch.size = min(batch.size, n.tot - cur.n), ...);
        rst     <- rbind(rst, cbind(BATCH = cur.b, Y = cur.rst));
        cur.n   <- cur.n + batch.size;
    }
    rst
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
    c(cur.test$estimate, cur.test$p.value, cur.test$conf.int, cur.test$conf.int[1] > p0);
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

sum.simon <- function(resp, design) {
    stop.1 <- sum(resp[1:design["n1"]]) <= design["r1"];
    rej    <- sum(resp) > design["r"];
    enroll <- ifelse(stop.1, design["n1"], design["n"]);
    mu     <- ifelse(stop.1, mean(resp[1:design["n1"]]), mean(resp));

    c(stop1  = stop.1,
      rej    = rej,
      enroll = enroll,
      rate   = mu)
}
