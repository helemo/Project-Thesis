library(SAVE)
#obtain dataframe
data("synthmodel", package = "SAVE")
data("synthfield", package = "SAVE")


#produce the object 
synth <- SAVE(response.name = "y", controllable.names = "x",
                calibration.names = "v", field.data = synthfield, ,
                model.data = synthmodel, mean.formula = ~ 1 + x,
                bestguess = list(v = 1.5))

#Obtain MCMC samples
set.seed(0)
synth <- bayesfit(object = synth,
                     prior = uniform(var.name = "v", lower = 0, upper = 3), n.iter = 20000)
plot(synth,option = "precision")
plot(synth,option = "calibration")

xnew <- data.frame(x = seq(from = 0.05, to = 3.05, length = 25))
xnew
valsynth <- validate(object = synth, newdesign = xnew, n.burnin = 100)

meansynth <- valsynth@validate[,"bias.corrected"]
upper_synth <- qnorm(0.90, valsynth@validate[,"bias.corrected"],valsynth@validate[,"tau.bc"])
lower_synth <- qnorm(0.10, valsynth@validate[,"bias.corrected"],valsynth@validate[,"tau.bc"])
length(lower_synth)
length(xnew[,1])
truemean <- 3.5*exp(-1.7*synthfield$x)+1.5
truemean

#Least squares estimate of calibration parameter
u_LSS = 0.63 
y_LSS <- 5*exp(-u_LSS*xnew)
y_LSSR <- 3.5*exp(-u_LSS*synthfield$x)+1.5
y_LSSR
#Computer model evaluated at posterior mean
u_postmean = 1.58
y_postmean <- 5*exp(-u_postmean*xnew)


plot(xnew[,1], meansynth,
     main="Predictions",
     ylab="y",
     xlab = "x",  
     ylim = c(0,5),
     type="l",
     col="black")
lines(xnew[,1],y_LSS[,1], col = "green")
#lines(synthfield$x,y_LSSR, col = "red")
lines(xnew[,1],y_postmean[,1], col = "blue")
lines(xnew[,1], lower_synth, col="black", lty=2)
lines(xnew[,1], upper_synth, col = "black", lty=2)
lines(synthfield$x, truemean, col="red")

points(synthfield$x, synthfield$y, col= "black",pch="+")
#legend("topright",
       #c("Est. mean w/90% quantiles", "Observed values","True mean"),
       #fill=c("green", "black", "blue")
#)

y_LSS[,1]
xnew[,1]
