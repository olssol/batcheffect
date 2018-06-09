rm(list=ls());
options(error = recover, warn=0);
source("batch_tools.R");

ARGS.INX <- as.numeric(commandArgs(trailingOnly=TRUE));
if (0 == length(ARGS.INX)) {
    ARGS.INX <- 1;
    print(paste("No argument entered"));
}

##-----SIMULATION SETTINGS---------

## hypothesis
P1        <- 0.30;
P0        <- 0.15;
ALPHA     <- 0.05; ##two sided
POWER     <- 0.8;

## response
MU        <- c(0, 2);
BETA      <- 0.3;

PAR.ERROR <- list(gamma   = list(error.type = "normal", ysig = 0.132),
                  delta   = list(error.type = "normal", ysig = 0, mu = 1),
                  epsilon = list(error.type = "normal", ysig = 1));

PAR.ERROR <- list(gamma   = list(error.type = "skewed", ysig = 0.132, skew.n = 2, skew.p = 0.5),
                  delta   = list(error.type = "normal", ysig = 0, mu = 1),
                  epsilon = list(error.type = "normal", ysig = 1));


test <- simu.error(n.simu = 10000, error.type = "skewed", ysig = 0.132, skew.n = 2,   skew.p = 0.5);
test <- simu.error(n.simu = 10000, error.type = "skewed", ysig = 0.132, skew.n = 1.5, skew.p = 0.1);
##plot(density(test));

## conduct
BATCH.SIZE <- 37;
N.BATCH    <- NULL;
## replications
NREPS      <- 2;

##-----SIMULATION ---------
rst <- simu.trial(p1 = P1, p0 = P0, batch.size = BATCH.SIZE, n.batch = N.BATCH,
                  mu = MU, beta = BETA,
                  par.error = PAR.ERROR,
                  alpha = ALPHA, power = POWER,
                  nreps = NREPS, seed = 10000);
