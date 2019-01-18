load("NEAM Q41 Jansen et al 2014 LGC model 10 km grid final model.RData")
##dfModelData

## Default inputs
## SPECIES  <- "Gadus morhua"
## QUARTER  <- 4
## KM       <- 50
## MINSIZE  <- 4
## MAXSIZE  <- 120
## MINYEAR  <- 2000
## MAXYEAR  <- 2015
## BY       <- 2
## DATFILE  <- "EBcod.RData"
## OUTFILE  <- paste0("results", QUARTER, ".RData")

## For scripting
## input <- parse(text=Sys.getenv("SCRIPT_INPUT"))
## print(input)
## eval(input)

## Load data
## library(DATRAS)
## d <- local({
##     load(DATFILE)
##     stopifnot(length(ls()) == 1)
##     get(ls())
## })
## stopifnot( class(d) == "DATRASraw" )

## Make grid
library(gridConstruct)
## grid <- gridConstruct(d,km=KM)
grid <- attr(dfModelData$position,"levelAttrib")$grid
## plot(grid)
## map("worldHires",add=TRUE)


## Data subset
d <- dfModelData
##d <- addSpectrum(d,cm.breaks=seq(MINSIZE,MAXSIZE,by=BY))
## d$haulid <- d$haul.id
## d <- subset(d, Quarter == QUARTER, Gear != "GRT")
## d <- subset(d, Year %in% MINYEAR:MAXYEAR )
## d <- subset(d, 25<HaulDur & HaulDur<35 )
## d <- as.data.frame(d)
d$sizeGroup <- factor(rep("0", nrow(d)))


library(mapdata)

## Set up time factor (careful with empty factor levels ! )
years <- as.numeric(levels(d$Year))
d$time <- factor(d$Year, levels=min(years):max(years))
d$haulid <- factor(1:nrow(d))
d$N <- d$NumbersRaised

## Set up spatial factor and grid
## grid <- gridConstruct(d,nearestObs=100)
##d$position <- gridFactor(d,grid)
Q0 <- -attr(grid,"pattern")
diag(Q0) <- 0; diag(Q0) <- -rowSums(Q0)
I <- .symDiagonal(nrow(Q0))

## Set up design matrix
## TODO: Gear by sizeGroup
##A <- sparse.model.matrix( ~ sizeGroup:time + Gear - 1, data=d)
##A0 <- sparse.model.matrix( ~ sizeGroup:time - 1, data=d)
##A0 <- sparse.model.matrix( ~ time - 1, data=d)
fixef <- ~ Year - 1 + I(log(GroundSpeed_kn)) + I(log(WingSpread_m/(1852))) + I(log(HaulDur_min/60))

## WARNING: old lgc package substitute missing values by zeros in design matrix !
## A <- sparse.model.matrix( fixef , data=d)
mf <- model.frame( ~ Year - 1 + I(log(GroundSpeed_kn)) + I(log(WingSpread_m/(1852))) + I(log(HaulDur_min/60)) , data=d, na.action=na.pass)
mm <- model.matrix(fixef,mf)
mm[is.na(mm)] <- 0
A <- as(mm,"dgCMatrix")

##A <- A[,-which.max(table(d$Gear)),drop=FALSE]

## NOTE: C++ assumes that the head(beta, nyears) is the year effect !!!
stopifnot(identical(colnames(A)[1:nlevels(d$Year)],
                    paste0("Year", levels(d$Year))))

##B <- cbind2(A,A0); B <- t(B)%*%B
B <- cbind2(A); B <- t(B)%*%B
if(min(eigen(B)$val)<1e-8)stop("Singular B")

data <- as.list(d)
data$A <- A
data$I <- I
data$Q0 <- Q0

data <- data[!sapply(data,is.character)]
data <- data[!sapply(data,is.logical)]
data$huge <- 0
data$h <- mean(summary(as.polygons(grid))$side.length)

stopifnot(nrow(data$A) == length(data$N))

library(TMB)
compile("model.cpp")
dyn.load(dynlib("model"))
obj <- MakeADFun(
    data=data,
    parameters=list(
        logdelta= 0 ,
        logkappa= 0 ,
        logkappa_static= 0 ,
        tphi_time= 0 ,
        logsigma= 0 ,
        beta= rep(0, ncol(A)) ,
        eta= array(0, c(nrow(Q0), nlevels(d$time) ) ),
        etanug= array(0, c(nlevels(d$haulid) ) ),
        eta_static= numeric(nrow(Q0))
        ##etamean = array(0, c(nlevels(d$sizeGroup) , nlevels(d$time)) )
        ),
    DLL="model",
    ##random=c("eta","etanug","beta")
    random=c("eta","etanug","eta_static","beta")
    )

print(obj$par)
runSymbolicAnalysis(obj)
system.time(obj$fn())

system.time(opt <- nlminb(obj$par, obj$fn, obj$gr))

system.time(sdr <- sdreport(obj))

lpb <- obj$env$last.par.best
pl <- obj$env$parList(par=lpb)
dim(pl$eta)
length(pl$eta_static)
I <- colMeans(sqrt(exp(t(t(pl$eta + pl$eta_static) + pl$beta[1:nlevels(d$Year)]))))

summary(sdr, "report")

est <- as.list(sdr, "Estimate")
std <- as.list(sdr, "Std. Error")
## Fixed effects
data.frame(Parameter=colnames(A), Estimate=est$beta, Std.Error=std$beta)

## Variance (Nugget)
exp(est$logsigma)^2

## Variance (YS)
Q <- Q0 + exp(est$logdelta) * I
exp(est$logkappa)^2 * mean(diag(TMB:::solveSubset(Q)))

## Variance (IS)
Q <- Q0 + exp(est$logdelta) * I
exp(est$logkappa_static)^2 * mean(diag(TMB:::solveSubset(Q)))
