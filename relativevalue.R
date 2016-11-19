library(chron)
library(zoo)

rm(list=ls(all=TRUE))

cdates <- function(maturity){

    startdate <- yearsbefore(maturity, 20)
    x <- seq.dates(startdate, maturity, by="months")
    cdates <- x[seq(1, length(x), by=6)] ; rm(x)

    return(cdates)
}

allbonddata <- function(dat.frame){

    bonddata2list <- function(dat.raw){

    maturity <- chron(dat.raw$Maturity, format = "y-m-d",
                      out.format = "m/d/y")
    yields <- zoo(dat.raw$Yield, today)
    coupon <- dat.raw$Coupon
    cdates <- cdates(maturity)
    prices <- zoo(bondprice(coupon, yields, cdates), today)

    return(list(ISIN = dat.raw$ISIN, maturity = maturity,
                issuer = as.character(dat.raw$Issuer),
                coupon = coupon, cdates = cdates,
                yields = yields, prices = prices))
    }

    bonddata <- vector("list", nrow(dat.frame))
    for(i in 1:nrow(dat.frame)){
        bonddata[[i]] <- bonddata2list(dat.frame[i,])
    }
    return(bonddata)
}

getparcurve <- function(issuer = NULL, type.i = NULL,  titl, today, spar.in = 0.65){
    if(!is.null(issuer)){
        bonddata.issuer <- bonddata.raw[bonddata.raw$Issuer %in% issuer,]
    }

    if(!is.null(type.i)){
        bonddata.issuer <- bonddata.raw[bonddata.raw$Type == type.i,]
    }

    maturities <- bonddata.issuer$Maturity
    maturities <- chron(maturities, format = "y-m-d", out.format = "m/d/y")
    bondrates <- bonddata.issuer$Yield
    issuers.set <- bonddata.issuer$Issuer
    issuers.set <-  zoo(issuers.set[!duplicated(maturities)],
                        maturities[!duplicated(maturities)])
    issuers.set <- issuers.set[index(issuers.set) > today + 366]
    bondrates <- zoo(bondrates[!duplicated(maturities)],
                     maturities[!duplicated(maturities)])
    bondrates <- bondrates[index(bondrates) > today + 366]
    # allrates <- c(cashrates, bondrates)
    allrates <- bondrates
    rate <- coredata(na.omit(allrates))
    tenor <- (index(na.omit(allrates)) - today) / 30
    plot(tenor, rate, main = titl,
         ylab = "Yield (per cent)", xlab = "Tenor (months)")
    ycurve <- smooth.spline(tenor, rate, spar = spar.in)
    ycurve.out <- predict(ycurve)$y
    lines(predict(ycurve, min(tenor):max(tenor)), lwd = 2.5, col = "blue")

    #shortlabels = c("CASH", "OIS1M", "OIS3M", "OIS6M")
    datalabels <- paste(as.character(coredata(issuers.set)),
                        as.character(index(allrates)))
    text(tenor, rate, datalabels, pos = 3, cex = .6,
         srt = 45, col = "darkgrey")
    relval <- round(100*(rate - ycurve.out), 2)
    results <- data.frame(bond=datalabels, yield=rate, 
                          model=round(ycurve.out,4), relval)
    return(results)
}



today <- chron("8/15/08")
# cashrates <- zoo(c(7.25, 7.095, 6.997, 6.922), today + c(1, 30, 90, 180))

klasses <- c("character", rep("factor", 2), "numeric")
bondlist <- read.csv("bondlist.csv", colClasses = klasses)
bondyields <- read.csv("bondyields.csv", colClasses = klasses[c(1, 4)])
bonddata.raw <- merge(bondlist, bondyields)

# ACGB relative value curve
getparcurve(today = today, spar.in = 0.5, type = "CGS", 
            titl = "ACGB yield curve")

# ACGB relative value curve, alternative method
getparcurve(today = today, spar.in = 0.5, issuer = "CTH", 
            titl = "ACGB yield curve")

# TCV curve
getparcurve(today = today, spar.in = 0.1, issuer= "TCV", 
            titl = "TCV yield curve")

# Top tier semi curve
getparcurve(today = today, spar.in = 0.65, issuer=c("QTC", "NSWTC", "TCV"), 
            titl = "TCORP, TCV, QTC yield curve")

# Lower tier semi curve
getparcurve(today = today, spar.in = 0.35, issuer=c("SAFA", "TASCORP", "WATC"), 
            titl = "SAFA, TASCORP, WATC yield curve")

# All supras (a bit crowded)
getparcurve(today = today, spar.in = 0.9, type = "Supra", 
            titl = "Supra yield curve")


