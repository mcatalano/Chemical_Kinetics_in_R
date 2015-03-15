## Set working directory and load package requirements
setwd("C:/Documents and Settings/Michael/Desktop/dA-AP EM/DNA Data/ODE Solve in R")

library(ggplot2)
library(reshape2)
library(deSolve)
library(minpack.lm)

## Read concentration/time data into data frame
df <- read.csv('KS35pH7.csv', header=F)
names(df) <- c('time','A','B','C','D')

## The rate function calculates results to the ODEs associated with the specified 
## kinetic model and returns the results as a list
rate <- function(t,c,params){

    ## Set rate constants to 'params' list elements passed to function
    k1 <- params$k1
    k2 <- params$k2
    k3 <- params$k3
    k4 <- params$k4
    k5 <- params$k5
    
    ## Create vector of calculated rates; 'c' is a vector of concentrations
    ## of each species at a given time point. c['x'] refers to an element of 'c'
    ## with name attribute 'x'
    r <- rep(0, length(c))
    r[1] <- -k1*c['A'] + k2*c['B'] - k3*c['A']
    r[2] <- -k2*c['B'] + k1*c['A']
    r[3] <- -k4*c['C'] + k3*c['A'] + k5*c['D']
    r[4] <- k4*c['C'] - k5*c['D']
    
    return(list(r))    
}

## The 'resid' function solves the ODEs specified by calling the 'rate' function in ode(), 
## then returns the residual between calculated and measured values
resid <- function(params){
    
    ## Set initial concentrations
    cinit <- c(A=100,B=0,C=0,D=0)
    
    ## Vector of time points
    t <- df$time
    
    ## k values set to 'params' elements passed to function
    k1 <- params[1]
    k2 <- params[2]
    k3 <- params[3]
    k4 <- params[4]
    k5 <- params[5]
    
    output <- ode(y=cinit, times=t, func=rate, parms=list(k1=k1,k2=k2,k3=k3,k4=k4,k5=k5))
    
    ## Store 'output' matrix as data frame
    outdf <- data.frame(output)
    
    ## Convert solution and measured data to long format
    preddf <- melt(outdf, id.var="time", variable.name="species", value.name="yield")
    expdf <- melt(df, id.var="time", variable.name="species", value.name="yield")
    
    ## Compute residual
    res <- preddf$yield - expdf$yield
    
    return(res)  
}

## Vector of estimated rate constants
params <- c(k1=.06, k2=.03, k3=.004, k4=.003, k5=.005)

## Least squares minimization of residual by adjusting specified parameters
fitval <- nls.lm(par=params, fn=resid)

## Store resulting parameters as a list
parest <- as.list(coef(fitval))

## Create a results matrix of solutions from least squares minimization by passing
## the solution parameters list, 'parest', to ode()

    ## Set initial concentrations
    cinit <- c(A=100,B=0,C=0,D=0)

    ## Create an arbitrarily incremented time sequence
    t <- seq(0, df[length(df$time),1], 0.1)
    
    ## Set parameters
    params <- as.list(parest)

    ## Solve
    output <- ode(y=cinit, times=t, func=rate, parms=params)

## Convert ODE solution matrix to labeled data frame
outdf <- data.frame(output)
names(outdf) <- c('time','A_pred','B_pred','C_pred','D_pred')

## Rearrange data to long format
pred <- melt(outdf, id.var=c('time'), variable.name='species', value.name='yield')
exp <- melt(df, id.var=c('time'), variable.name='species', value.name='yield')

## Open graphics device
png(file = "Time_course.png")

## Plot results
p <- ggplot(data=pred, aes(x=time, y=yield, color=species)) + geom_line()
p <- p + geom_point(data=exp, aes(x=time, y=yield, color=species))
p <- p + labs(x='Time (hr)', y='% Yield')

print(p)

## Close device
dev.off()