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
P1        <- 0.50;
P0        <- 0.15;
ALPHA     <- 0.05; ##two sided
POWER     <- 0.8;

## response
MU        <- 0;
PAR.ERROR <- list(gamma   = list(error.type = "normal", ysig = 0.1),
                  delta   = list(error.type = "normal", ysig = 0.1, mu = 1),
                  epsilon = list(error.type = "normal", ysig = 1));

## conduct
BATCH.SIZE <- 9;
N.BATCH    <- NULL;
## replications
NREPS <- 1000;


##-----SIMULATION ---------
rst <- simu.trial(p1 = P1, p0 = P0, batch.size = BATCH.SIZE, n.batch = N.BATCH,
                  mu = MU, par.error = PAR.ERROR,
                  alpha = ALPHA, power = POWER,
                  nreps = NREPS, seed = 10000);
rst

##-----generate table ---------
sub.txt.w.no(1:100, template.f = "temp.txt", out.f = "rst.txt");
