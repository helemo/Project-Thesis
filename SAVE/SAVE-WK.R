library(ggplot2)
library(reticulate)
library(SAVE)
library(DiceKriging)
library(lhs)
#library(laGP)

np = import("numpy")

#setwd("Z:/PROSJEKTOPPGAVE/Script")
#setwd("C:/Users/hming/OneDrive - NTNU/Skole/5.klasse/PROSJEKTOPPGAVE")


### Michail WK3 Model ###

# Wk2 pressure simualtor
pressurecalculator_wk2 = function(flow,time,R,C){
  library(reticulate)
  np = import("numpy")
  len = function(x) length(x)
  
  impedance = function(w,R,C){# Computes the electrical impedance Z(omega) for the 2 element Windkessel
    Z=R/(1.0+1i*w*R*C) #w = 2pif?
    return(Z)
  }
  real_irfft = function(Z,w,t){
    N=len(t)
    Nf=len(w)
    Z=Z/N
    
    r=np$zeros(len(t)) + np$real(Z[1])
    for (k in 1:(Nf-1)){
      r = r + np$real(2*Z[k+1]*np$exp(1i*w[k+1]*t))
    }
    return(r)
  }
  
  TT = time[len(time)] - time[1]
  
  #Load data
  Q = np$asarray(flow)
  t = time
  dt = t[2:len(t)] - t[1:(len(t)-1)]
  #Interpolate if uneven spacing
  if (np$max(dt) > 1.1*np$min(dt)){
    t_new = np$linspace(0.0, time[-1], len(time))
    q_new = np$interp(t_new,time,Q)
    
    time = t_new
    Q = q_new
  }
  
  Qfft = np$fft$rfft(Q) #compute 1D Fourier Trans for Real input 
  
  #w = np$fft$rfft(t_mark)*2.0*np$pi
  
  w = 2.0*np$pi*np$fft$rfftfreq(len(time),time[2]-time[1]) # return the discrete FT sample freq = 2pif!!!
  P = real_irfft(impedance(w,R,C)*Qfft,w,t)
  
  #P = np$fft$irfft(modref.impedance(w)*Qfft,n=len(t))
  #print("Lengde: ", len(t_mark), len(w), len(Qfft), len(P))
  
  P_sys = np$max(P)
  P_dia = np$min(P)
  PP = P_sys-P_dia
  SV = np$trapz(Q,x=time) #integrate along the given axis using the composite trapezoidal rule
  CO = SV/T*60.0/1000.0 #L/min
  P_br_sys = P_sys/0.9 #3.0*(P_sys+P_dia)/2.0 - 2.0*P_dia
  MAP = CO/60.0*1000.0*R
  SW = MAP*SV
  
  l = list(P=P, time = time, P_sys = P_sys, P_dia = P_dia)
  return(l)
}

# WK3 pressure simulator
pressurecalculator_wk3 = function(flow,time,R,C,Z){
  library(reticulate)
  np = import("numpy")
  len = function(x) length(x)
  
  impedance = function(w,R,C,Z_ao){
    Z=(Z_ao + R + 1i*w*R*C*Z_ao)/(1.0+1i*w*R*C)
    return(Z)
  }
  
  real_irfft = function(Z,w,t){
    N=len(t)
    Nf=len(w)
    Z=Z/N
    
    r=np$zeros(len(t)) + np$real(Z[1])
    for (k in 1:(Nf-1)){
      r = r + np$real(2*Z[k+1]*np$exp(1i*w[k+1]*t)) # check indices again
    }
    return(r)
  }
  
  TT = time[len(time)] - time[1] #range of time-values
  Q = np$asarray(flow) # Load data
  
  t = time
  dt = t[2:len(t)] - t[1:(len(t)-1)]
  
  #Interpolate if uneven spacing
  if (np$max(dt) > 1.1*np$min(dt)){
    t_new = np$linspace(0.0, time[-1], len(time))
    q_new = np$interp(t_new,time,Q)
    
    time = t_new
    Q = q_new
  }
  
  
  Qfft = np$fft$rfft(Q) #Fourier transform of real part 
  #w = np$fft$rfft(t_mark)*2.0*np$pi
  w = 2.0*np$pi*np$fft$rfftfreq(len(time), time[2]-time[1])
  
  #own check!
  Z_est = impedance(w,R,C,Z)
  
  P = real_irfft(impedance(w,R,C,Z)*Qfft,w,t) 
  #P = np$fft$irfft(modref$impedance(w)*Qfft,n=len(t))
  #print("Lengde: ", len(t_mark), len(w), len(Qfft), len(P))
  P_sys = np$max(P)
  P_dia = np$min(P)
  PP = P_sys-P_dia
  SV = np$trapz(Q,x=time) #integrate along the given axis using the composite trapezoidal rule
  CO = SV/TT*60.0/1000.0 #L/min
  P_br_sys = P_sys/0.9 #3.0*(P_sys+P_dia)/2.0 - 2.0*P_dia
  MAP = CO/60.0*1000.0*R
  SW = MAP*SV
  
  l = list(P=P, time = time, P_sys = P_sys, P_dia = P_dia, Z_est=Z_est)
  return(l)
}

create_field_wk3 = function(flow, time, Rreal, Creal, Zreal, sdnoise, cycles, ind, seed=123){
  Pwk3 = pressurecalculator_wk3(flow,time,Rreal,Creal,Zreal)$P
  l = cycles*length(ind) #run nr of cycles per chosen indice
  df = data.frame(y = rep(NA, l), t = rep(time[ind],cycles), Q = rep(flow[ind], cycles))
  #set.seed(seed)
  y = matrix(NA, nrow = length(ind), ncol = cycles)
  for(i in 1:cycles)  y[,i] = Pwk3[ind] + rnorm(length(ind), 0, sdnoise)
  df$y = as.vector(y)
  return(list(field_data = df, Pwk3 = data.frame(t = time,P = Pwk3)))
}


###########################################################################

Flow_dat <- readRDS("WK2DatNew.rds")
flow <- Flow_dat[,2]*0.7

#Standardise flow
mean <- mean(flow)
sd <- sd(flow)
flow_std <- (flow-mean)/sd

time <- Flow_dat[,1]
#plot(time,flow,xlab = "time(s)",ylab="flow(ml/min)")

#Run for HP and nonHP
R <- 1.4
C <- 1.3
Z_HP <- 0.1
Z_nonHP <- 0.02

cycles <- 10 
sdnoise2 <- 4

#want to sample unevenly in time due to the larger curvature
ind_model <-c(round(seq(1,12,length.out = 4)),round(seq(20,50,length.out = 6)),round(seq(55,101,length.out = 6)))


#Build emulator 
#Create sample grid of R and C
set.seed(1234)
RC <- maximinLHS(cycles,2)*2.5 + 0.5

flow_sim <- flow[ind_model]
time_sim <- time[ind_model]

pM <- c()
RM <- c()
CM <- c()
QM <- c()
tM <- c()

for (i in 1:cycles){
  RM <- c(RM, rep(RC[i,1],length(ind_model)))
  CM <- c(CM, rep(RC[i,2],length(ind_model)))
  QM <- c(QM, flow_sim)
  tM <- c(tM, time_sim)
  
  WK_cycle <- pressurecalculator_wk2(flow,time, RC[i,1],RC[i,2])
  
  pM <- c(pM,WK_cycle$P[ind_model])
  
}

df_pM <- data.frame(P = pM)

#standardize Q for GP fit
mu_Q <- mean(QM)
sd_Q <- sd(QM)
stndrd_Q <- (QM-mu_Q)/sd_Q

df_xM <- data.frame(R = RM, C = CM, Q = stndrd_Q, t = tM)
mat_xM <- as.matrix(df_xM)

#define surrogate fit size for SAVE
ind_SAVE <- round(seq(1,160,length.out = 80))

#BUILD EMULATOR FOR CHECK
#Create a diceKriging object, estimating the unknown parameters using MLE 
test_model <- km(design = mat_xM[ind_SAVE,], response = df_pM[ind_SAVE,], nugget = 1e-8, multistart = 100,covtype = "powexp" )#,coef.var = sdnoise2)
test_model@covariance@sd2

#parameter estimate variaions?
test_model@covariance@range.val
png("ParVarMLE-4.png")
plot(test_model@covariance@range.val[1:2])
text(c(1.2,2),test_model@covariance@range.val[1:2],as.character(test_model@covariance@range.val[1:2]))
dev.off()

#parameter estimates
test_model@covariance@shape.val
png("ParEstMLE-4.png")
plot(test_model@covariance@shape.val[1:2])
text(c(1.2,2),test_model@covariance@shape.val[1:2],as.character(test_model@covariance@shape.val[1:2]))
dev.off()

#show(test_model)

#all indices defined
#new_sample <- round(seq(1,101, length.out = 101))
flow_new <- flow_std
flow_new
time_new <- time
time_new
df_new_X <- data.frame(R = rep((R),length(flow_new)), C = rep(C, length(flow_new)), Q = flow_new, t = time_new)

y_hat <- predict(test_model, as.matrix(df_new_X),"SK")

var <- y_hat$sd^2
cat("Var P: \ ", var," \ ")
var
#cat('variance of predicted values',var)

#16 indices
#WK_compare <- pressurecalculator_wk2(flow, time, R, C)$P[ind_model]
#all indices 
WK_compare <- pressurecalculator_wk2(flow, time, R, C)$P

#NB remember to change this x = time/time_sim when you change indices!!!!!!
df_plot <- data.frame(x = time_new, y1 = y_hat$mean,y2 = WK_compare, lb = y_hat$lower95, ub =y_hat$upper95)

ggplot(df_plot,aes(x = Time(s),y = PredictedPressure))+
  geom_line(aes(x = x, y = y1, color = "My Surrogate Fit"))+
  geom_line(aes(x = x, y = y2, color = "WK2 model"))+
  geom_line(aes(x = x, y = ub, color = "upper 95"), linetype = "dashed")+
  geom_line(aes(x = x, y = lb, color = "lower 95"), linetype = "dashed")+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  ggsave("SurrogateFit2-4.png")

# SAVE fit - same model design for both to ensure numerical stability
#use setup that makes a good emulator from emulator test
SAVE_ind_model <- ind_SAVE #c(round(seq(1,150,length.out = 8)),round(seq(160,320,length.out = 8)))#to be 16 in total but sample from whole space
#use same as for model
SAVE_ind_field <- round(seq(1,160,length.out = 16)) #field is form 101 field obs - need 16 field x 20 cycles

#create model obs
Model <- cbind(mat_xM, pM)
colnames(Model)[5] <- "P"
SAVE_model <- Model[SAVE_ind_model,]

#create field obs 
WK_field <- create_field_wk3(flow,time,R,C,Z_nonHP, cycles = cycles, ind = ind_model, sdnoise = 4) #ind is the same amount of field obs obtained in emulator fit
plot(WK_field$field_data$t , WK_field$field_data$y)
max(WK_field$field_data$y)

#standardize Q for field
mu_Q_field <- mean(WK_field$field_data$Q)
sd_Q_field <- sd(WK_field$field_data$Q)
stndrd_Q_field <- (WK_field$field_data$Q-mu_Q_field)/sd_Q_field

#create df of field data
df_field <- data.frame(Q = stndrd_Q_field, t = WK_field$field_data$t,P = WK_field$field_data$y)

Field <- as.matrix(df_field)
SAVE_field <- Field[SAVE_ind_field,]

kc = list(multistart = 10^2) # multistart optimization for the DiceKrig model
tic1 = Sys.time()
set.seed(123)

#create save object - emulator
pw <- SAVE(response.name = "P",
           controllable.names = c("Q","t"),
           calibration.names = c("R","C"), field.data = SAVE_field,
           model.data = SAVE_model, mean.formula = ~ 1, bestguess = list(R = 1.4,C=1.3),kriging.controls = kc)


#summary(pw)

# fit emulator - relating reality to computer model
set.seed(1444)
pwbayes <- bayesfit(object = pw, prior = c(uniform("R", upper = 3,lower = 0.5),uniform("C", upper = 3,lower = 0.5)), n.iter = 1e7, n.burnin = 10000, n.thin = 2000) #prev 20 000, 100,2

png("Calibration-4.png")
plot(pwbayes, option = "calibration")
dev.off()
plot(pwbayes, option = "calibration")


png("Precision-4.png")
plot(pwbayes, option = "precision")
dev.off()
plot(pwbayes, option = "precision")



# predict with calibration params at the posterior mean
#All indices
#xnew <- data.frame(Q = flow,t = time)
#selfdefined indices
xnew <- data.frame(Q = flow_new,t = time_new)
xnew
set.seed(1444)
valpw <- validate(object = pwbayes, newdesign = xnew, calibration.value = "mean", n.burnin = 1000)

plot(valpw)

png("Validation-4.png")
plot(valpw)
dev.off()

#calculate 95% pred int of the book example using save 
upper_SAVE <- qnorm(0.95, valpw@validate[,"bias.corrected"],valpw@validate[,"tau.bc"])
lower_SAVE <- qnorm(0.05, valpw@validate[,"bias.corrected"],valpw@validate[,"tau.bc"])

mean <- valpw@validate[,"bias.corrected"]

#pure model estimate for comparison
prpw <- predictreality(pwbayes,xnew)
model <- prpw@modelpred
mean_pred <- apply(model,2,mean)

#Change this from time/time[SAVE_ind_field] amd wKfield[SAVE_ind_field] ind when using different indices!!!!
df_SAVE <- data.frame(x = time_new, y1 = mean, lb = lower_SAVE, ub = upper_SAVE, true = WK_field$Pwk3$P, pure_model = mean_pred, comp = WK_compare)

MSE <- sum((WK_field$Pwk3$P-mean)^2)/101
png("MSE-4.png")
plot(MSE)
dev.off()

ggplot(df_SAVE, aes(x = Time(s), y = AorticPressure(mmHg)))+
  geom_point(data = df_field, aes(x = t, y = P),color = "grey")+
  geom_line(aes(x = x, y = y1, color = "WK2 Bias Corrected"))+
  geom_line(aes(x = x, y = true, color = "True Mean"))+
  geom_line(aes(x = x, y = comp, color = "WK2 True Parameter"))+
  geom_line(aes(x = x, y = lb), color = "#33CCCC", size = 0.3, linetype = "dashed") + 
  geom_line(aes(x = x, y = ub), color = "#33CCCC", size = 0.3, linetype = "dashed")+
  geom_line(aes(x = x, y = pure_model, color = "Pure Model"))+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  ggsave("SAVEfit-4.png")

