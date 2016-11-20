NelsonSiegel =
function(Yield, Maturity)
{    # A function written by Diethelm Wuertz
   
    # Description:
    #    Fit the Yield Curve by the Nelson-Siegel Method
    #
    # Details:
    #    This function finds a global solution. The start values for the
    #    betas are solved exactly as a function of tau using OLS.
    #
    # Copyright:
    #    Diethelm Wuertz, (c) 2004 fBonds
    #
    # Source:
    #    Partial copy from 'fBonds' from 'Rmetrics' (unpublished).
   
    # FUNCTION:
   
    # Find Optimal Start Solution by OLS of beta's vs. Yields:
    n = length(Maturity)
    gmin = 1.0e99
    for (i in 1:n) {    
        tau = Maturity[i]
        x = Maturity/tau        
        a = matrix(rep(NA, times = 9), nrow = 3)
            a[1,1] = 1
            a[1,2] = a[2,1] = mean((1-exp(-x))/x)
            a[1,3] = a[3,1] = mean((1-exp(-x))/x - exp(-x))
            a[2,2] = mean( ((1-exp(-x))/x)^2 )
            a[2,3] = a[3,2] = mean(((1-exp(-x))/x)*((1-exp(-x))/x-exp(-x)))
            a[3,3] = mean(((1-exp(-x))/x - exp(-x))^2)
        b = c(
            mean ( Yield ),
            mean ( Yield *  ((1-exp(-x))/x)),
            mean ( Yield * (((1-exp(-x))/x - exp(-x)))))            
        beta = solve(a, b)
        yfit = beta[1] + beta[2]*exp(-x) + beta[3]*x*exp(-x)
        fmin = sum( (Yield-yfit)^2 )
        if (fmin < gmin) {
            gmin = fmin
            gvec = c(beta, tau)
        }
    }
                   
    # Function to be optimized:       
    fx <- function(Maturity, x) {
        x[1] + x[2] * (1-exp(-Maturity/x[4]))/(Maturity/x[4]) +
        x[3] * 
((1-exp(-Maturity/x[4]))/(Maturity/x[4])-exp(-Maturity/x[4]))
    }
    func <- function(x) { sum( (Yield - fx(Maturity, x))^2 ) }

    # Optimize:
    fit = nlminb(objective = func, start = gvec)
    fit$start = gvec
    names(fit$par) = c("beta1", "beta2", "beta3", "tau")
       
    # Plot Curve:
    yfit = fx(Maturity, gvec)
    plot(Maturity, Yield, ylim = c(min(c(Yield, yfit)), max(c(Yield, 
yfit))),
        pch = 19, cex = 0.5, main = "Nelson-Siegel" )
    lines(Maturity, yfit, col = "steelblue")
           
    # Return Value:
    fit
}



NelsonSiegel2 =
function(Yield, Maturity, tau)
{    # A function written by Diethelm Wuertz
   
    # Description:
    #    Fit the Yield Curve by the Nelson-Siegel Method
    #
    # Details:
    #    This function finds a global solution. The start values for the
    #    betas are solved exactly as a function of tau using OLS.
    #
    # Copyright:
    #    Diethelm Wuertz, (c) 2004 fBonds
    #
    # Source:
    #    Partial copy from 'fBonds' from 'Rmetrics' (unpublished).
   
    # FUNCTION:
   
    # Find Optimal Start Solution by OLS of beta's vs. Yields:
    n = length(Maturity)
    gmin = 1.0e99
    for (i in 1:n) {    
        tau = Maturity[i]
        x = Maturity/tau        
        a = matrix(rep(NA, times = 9), nrow = 3)
            a[1,1] = 1
            a[1,2] = a[2,1] = mean((1-exp(-x))/x)
            a[1,3] = a[3,1] = mean((1-exp(-x))/x - exp(-x))
            a[2,2] = mean( ((1-exp(-x))/x)^2 )
            a[2,3] = a[3,2] = mean(((1-exp(-x))/x)*((1-exp(-x))/x-exp(-x)))
            a[3,3] = mean(((1-exp(-x))/x - exp(-x))^2)
        b = c(
            mean ( Yield ),
            mean ( Yield *  ((1-exp(-x))/x)),
            mean ( Yield * (((1-exp(-x))/x - exp(-x)))))            
        beta = solve(a, b)
        yfit = beta[1] + beta[2]*exp(-x) + beta[3]*x*exp(-x)
        fmin = sum( (Yield-yfit)^2 )
        if (fmin < gmin) {
            gmin = fmin
            gvec = c(beta, tau)
        }
    }
                   
    # Function to be optimized:       
    fx <- function(Maturity, x) {
        x[1] + x[2] * (1-exp(-Maturity/tau))/(Maturity/tau) +
        x[3] * 
((1-exp(-Maturity/tau))/(Maturity/tau)-exp(-Maturity/tau))
    }
    func <- function(x) { sum( (Yield - fx(Maturity, x))^2 ) }

    # Optimize:
    fit = nlminb(objective = func, start = gvec)
    fit$start = gvec
    names(fit$par) = c("beta1", "beta2", "beta3", "tau")
       
    # Plot Curve:
    yfit = fx(Maturity, gvec)
    plot(Maturity, Yield, ylim = c(min(c(Yield, yfit)), max(c(Yield, 
yfit))),
        pch = 19, cex = 0.5, main = "Nelson-Siegel" )
    lines(Maturity, yfit, col = "steelblue")
           
    # Return Value:
    fit
}


nsplot <- function(rates, today = chron("8/12/08"), excl = NULL){
  startcurve <- cbind(coredata(rates), (index(rates) - today) / 30)
  if(!is.null(excl)){
    NelsonSiegel2(startcurve[-excl,1], startcurve[-excl,2], tau = 10)
  }
  else{
    NelsonSiegel2(startcurve[,1], startcurve[,2], tau = 10)
  }
}

# nsplot(zeroes$startval)
# nsplot(zeroes$startval, excl = c(1, 2, 4))
# nsplot(zeroes$zeroes.annual)





