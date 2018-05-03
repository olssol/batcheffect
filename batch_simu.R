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
P1        <- 0.3;
P0        <- 0.15;
ALPHA     <- 0.05; ##two sided
POWER     <- 0.8;

## response
MU        <- 0;
PAR.ERROR <- list(gamma   = list(error.type = "normal", ysig = 1),
                  delta   = list(error.type = "normal", ysig = 1, mu = 1),
                  epsilon = list(error.type = "normal", ysig = 1));
truth.pts <- simu.all.pts(100000, batch.size = 1, mu = MU, par.error = PAR.ERROR)[,"Y"];
CUT.Y.P1  <- quantile(truth.pts, 1-P1);
CUT.Y.P0  <- quantile(truth.pts, 1-P0);

## design
DESIGN.SIMON  <- get.simon.2.opt(pu = P0, pa = P1, ep1 = ALPHA,   ep2 = 1 - POWER);
DESIGN.SINGLE <- get.single.size(p0 = P0, p1 = P1, alpha = ALPHA, beta = 1 - POWER);

## conduct
BATCH.SIZE <- 3;

## replications
NREPS <- 1000;

##-----SIMULATE---------
all.rst <- NULL;
for (rep in 1:NREPS) {
    print(rep);

    ## ---- simon design ------------
    cur.simon.pt <- simu.all.pts(n.tot      = DESIGN.SIMON["n"],
                                 batch.size = BATCH.SIZE,
                                 mu         = MU,
                                 par.error  = PAR.ERROR);
    simon.h0 <- sum.simon(cur.simon.pt[,"Y"] > CUT.Y.P0, DESIGN.SIMON);
    simon.h1 <- sum.simon(cur.simon.pt[,"Y"] > CUT.Y.P1, DESIGN.SIMON);

    ##----single stage design -----
    cur.single.pt <- simu.all.pts(n.tot     = DESIGN.SINGLE["n"],
                                 batch.size = BATCH.SIZE,
                                 mu         = MU,
                                 par.error  = PAR.ERROR);
    single.h0 <- sum.single(cur.single.pt[,"Y"] > CUT.Y.P0, ALPHA, P0);
    single.h1 <- sum.single(cur.single.pt[,"Y"] > CUT.Y.P1, ALPHA, P0);

    all.rst   <- rbind(all.rst, c(simon.h0, simon.h1, single.h0, single.h1));
}
