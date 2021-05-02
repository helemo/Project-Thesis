library(SAVE)
library(GGally)
setwd("C:/Users/hming/OneDrive - NTNU/Skole/5.klasse/PROSJEKTOPPGAVE/Examples")

ball <- read.csv("ball.csv")
ggpairs(ball)


#predict time for ball to hit ground for different heights
#fit an isotropic GP directly to the data y_hat_F
library(laGP)
field.fit <- newGP(as.matrix(ball$height), ball$time, d=0.1, 
                   g=var(ball$time)/10, dK=TRUE)
eps <- sqrt(.Machine$double.eps)

#perform joint inference by iterating over marginals of lenghtscale and nugget
mle <- jmleGP(field.fit, drange=c(eps, 10), grange=c(eps, var(ball$time)),
              dab=c(3/2, 8))

#USING FIRST MODEL TO TEST THIS PACKAGE
#consider predictions on a testing grid wobbleX
#hs grid of heights in terms of coded inputs
hr <- range(ball$height)
hs <- seq(0, 1, length=100)
heights <- hs*diff(hr) + hr[1]
p <- predGP(field.fit, as.matrix(heights), lite=TRUE)
p$mean
deleteGP(field.fit)

par(mfrow=c(1,1))
#plot GP fit to field data, red dashed is predictive sd
plot(ball, xlab="height", ylab="time")
lines(heights, p$mean, col=4)
lines(heights, qnorm(0.05, p$mean, sqrt(p$s2)), lty=2, col=4)
lines(heights, qnorm(0.95, p$mean, sqrt(p$s2)), lty=2, col=4)
lines(heights, 10*sqrt(p$s2)-0.6, col=2, lty=3, lwd=2)
legend("topleft", c("Fhat summary", "Fhat sd"), lty=c(1,3), 
       col=c(4,2), lwd=1:2)

#coupling with known physical model t = sqrt(2h/g)
#g is calibration param - our u
#mapping natural inputs [h,g] to coded [x,u] in [0,1]^2
timedrop <- function(x, u, hr, gr) 
{
  g <- diff(gr)*u + gr[1]
  h <- diff(hr)*x + hr[1]
  return(sqrt(2*h/g))
}

#aquiring a prior for p(u) for g 
gr <- c(6, 14)

#21 unique heights, construct LHS in 2d and perform simulations at those locations
library(lhs) 
XU <- maximinLHS(21, 2) 
XU
yM <- timedrop(XU[,1], XU[,2], hr, gr)

#renaming and converting back to natural grid
time <- yM
height <- XU[,1]*diff(hr)+hr[1]
u <- XU[,2]*diff(gr)+gr[1]

df <- data.frame(height,time,u)
length(df)
  
#obtain dataframe
sample_field <- as.integer(seq(1,60,length.out = 15))
sample_model <- as.integer(seq(1,21,length.out = 20))

ballfield <- read.csv("ball.csv")[sample_field,]
plot(ballfield)
#sample evenly spaced
#ballfield = ballfield[seq(1,60,length.out = 10),]
ballmodel = df[1:10,]
ballmodel2 = df[sample_model,]
#create object using SAVE 
bs <- SAVE(response.name = "time", controllable.names = "height",
                calibration.names = "u", field.data = ballfield, 
                model.data = ballmodel2, mean.formula = ~ 1,
                bestguess = list(u = 10))

summary(bs)

#Obtain MCMC samples
set.seed(0)
bs <- bayesfit(object = bs,
                  prior = uniform(var.name = "u",upper=14,lower=6), n.iter = 200,n.burnin = 100, n.thin = 2,verbose=TRUE)
plot(bs,option = "precision")
plot(bs,option = "trace")
plot(bs,option = "calibration")
#length(bs@mcmcsample)
bs@bestguess




#predicitons of reality, assess pure-model and bias corrected
heightnew <- seq(from = 0, to = 4.5,  length=80)
newnew <- expand.grid(height = heightnew)
#xnew <- expand.grid(current = curr, load = load, thickness = g)
set.seed(0)
valbs <- validate(object = bs, newdesign = newnew,
                  calibration.value = "mean", n.burnin = 2)



(plot(valbs))
summary(valbs)

par(mfrow=c(1,1))

#calculate 90% pred int of the book example using save 
upper_SAVE <- qnorm(0.90, valbs@validate[,"bias.corrected"],valbs@validate[,"tau.bc"])
lower_SAVE <- qnorm(0.10, valbs@validate[,"bias.corrected"],valbs@validate[,"tau.bc"])

mean <- valbs@validate[,"bias.corrected"]

#fetch book example results only using KOH bayesian estimation
mean_book <- pm$mean + pb$mean
Sigma.Rhat <- pm$Sigma + pb$Sigma - diag(eps, length(m))
s_book <- sqrt(diag(Sigma.Rhat))

upper_book <- qnorm(0.9, m, s)
lower_book <- qnorm(0.1, m, s)

#The simple fnc for calculation of time for ball to hit the ground
analytical <- sqrt(2*heightnew/9.81)
upper_analytical <- qnorm(0.90,analytical,0.1)
lower_analytical <- qnorm(0.1,analytical,0.1)

length(heightnew)
#sd_SAVE <- valbs@validate[,"tau.bc"]
#BURDE JEG PLOTTE SD SOM DE GJORDE  IBOKA

par(mfrow=c(1,1))
plot(heightnew, mean,
     main="Prediction mean ",
     ylab="Time",
     xlab = "Height",  
     ylim = c(0.2,1.2),
     type="l",
     col="blue")
lines(heightnew, lower_SAVE, col="blue", lty=2)
lines(heightnew, upper_SAVE, col = "blue", lty=2)
#lines(heightnew, analytical, col = "green")
lines(heights, mean_book, col="red")
lines(heights, lower_book, col="red",lty=2)
lines(heights, upper_book, col = "red", lty=2)
#lines(heightnew, lower_analytical, col="red", lty=2)
#lines(heightnew, upper_analytical, col = "red",lty=2)
points(ball$height, ball$time, col= "black")
#legend("topleft",
       #c("SAVE", "Field","Book"),
       #fill=c("blue","black","green")
#)



############################################

