library(ggplot2)
library(reticulate)
library(SAVE)
library(DiceKriging)
library(lhs)
#library(laGP)

np = import("numpy")

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


Flow_dat <- readRDS("WK2DatNew.rds")
flow <- Flow_dat[,2]*0.7
mean <- mean(flow)
sd <- sd(flow)

#standardize flow and check if new mean = 0 and sd = 1
flow_std <- (flow-mean)/sd
time <- Flow_dat[,1]

#plot(time,flow_std,xlab = "time(s)",ylab="flow(ml/min)")

#Run for HP and nonHP
R <- 1.4
C <- 1.3
Z_HP <- 0.1
Z_nonHP <- 0.02

cycles <- 20 #20 for sim!
sdnoise2 <- 4
ind_model <-c(round(seq(1,12,length.out = 4)),round(seq(20,50,length.out = 6)),round(seq(55,101,length.out = 6)))
ind_model
#Build emulator 

#Create sample grid of R and C
set.seed(1234)
RC <- maximinLHS(cycles,2)*2.5 + 0.5
RC

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

#subtract mean from P
mu_p <- mean(pM)
#df_pM <- data.frame(P = (pM-mu_p))

#use regular pressure
df_pM <- data.frame(P = pM)

#standardize Q 
mu_Q <- mean(QM)
sd_Q <- sd(QM)

stndrd_Q <- (QM-mu_Q)/sd_Q

#create design matrix
df_xM <- data.frame(R = RM, C = CM, Q = stndrd_Q, t = tM)
mat_xM <- as.matrix(df_xM)
#mat_xM

#ind_SAVE <- round(seq(1,160,length.out = 80))
'ind_SAVE <- c(round(seq(1,5,length.out = 2)),
              round(seq(1,5,length.out = 1))+16*7,
              round(seq(1,5,length.out = 2))+16*5,
              round(seq(1,5,length.out = 1))+16*9,
              round(seq(1,5,length.out = 2))+16*2,
              round(seq(1,5,length.out = 1))+16*6,
              round(seq(1,5,length.out = 2))+16*3,
              round(seq(1,5,length.out = 1))+16*8,
              round(seq(1,5,length.out = 2))+16*4,
              
              
              round(seq(8,12,length.out = 2))+16*7,
              round(seq(8,12,length.out = 1))+16*5,
              round(seq(8,12,length.out = 2))+16*9,
              round(seq(8,12,length.out = 1))+16*2,
              round(seq(8,12,length.out = 2))+16*6,
              round(seq(8,12,length.out = 1))+16*3,
              round(seq(8,12,length.out = 2))+16*8,
              round(seq(8,12,length.out = 1))+16*4,
              
              
              round(seq(12,16,length.out=2)),
              round(seq(12,16,length.out=1))+16*7,
              round(seq(12,16,length.out=2))+16*5,
              round(seq(12,16,length.out=1))+16*9,
              round(seq(12,16,length.out=2))+16*2,
              round(seq(12,16,length.out=1))+16*8,
              round(seq(12,16,length.out=2))+16*4
              
              
              )'
ind_SAVE <- round(seq(1,160,length.out = 320))
ind_SAVE
length(ind_SAVE)




#BUILD EMULATOR FOR CHECK
#Create a diceKriging object, estimating the unknown parameters using MLE 
test_model <- km(design = mat_xM[ind_SAVE,], response = df_pM[ind_SAVE,], nugget = 1e-8, multistart = 100,covtype = "powexp" )#,coef.var = sdnoise2)

test_model@covariance@sd2

#parameter estimate variaions?
test_model@covariance@range.val
png("ParVarMLE.png")
plot(test_model@covariance@range.val[1:2])
text(c(1.2,2),test_model@covariance@range.val[1:2],as.character(test_model@covariance@range.val[1:2]))
dev.off()

#parameter estimates
test_model@covariance@shape.val
png("ParEstMLE.png")
plot(test_model@covariance@shape.val[1:2])
text(c(1.2,2),test_model@covariance@shape.val[1:2],as.character(test_model@covariance@shape.val[1:2]))
dev.off()

show(test_model)

y_sim <- simulate(test_model)
yc <- c(y_sim)
yc
df_chck <- data.frame(x = mat_xM[ind_SAVE,4], y1 = yc)

ggplot(df_chck,aes(x = Time(s),y = Realisations))+
  geom_point(aes(x = x, y = y1))+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  ggsave("Realisations.png")

#estimated R and C values  
#df_new_X <- data.frame(R = rep(R,length(flow_sim)), C = rep(C, length(flow_sim)), Q = flow_sim, t = time_sim)
#all indices
flow_new <- flow_std
flow_new
time_new <- time
time_new
df_new_X <- data.frame(R = rep(R,length(flow_new)), C = rep(C, length(flow_new)), Q = flow_new, t = time_new)

y_hat <- predict(test_model, as.matrix(df_new_X),"SK")

y_hat_mean <- mean(y_hat$mean)
y_hat_sd <- y_hat$sd
test_sd <- sd(y_hat$mean)


#16 indices
#WK_compare <- pressurecalculator_wk2(flow, time, R, C)$P[ind_model]
#all indices 
WK_compare_optimal <- pressurecalculator_wk2(flow, time, R, C)$P
mu_WKcomp <- mean(WK_compare_optimal)

#test for a range of RC values
WK_compare_1 <- pressurecalculator_wk2(flow, time, RC[1,1], RC[1,2])$P
WK_compare_5 <- pressurecalculator_wk2(flow, time, RC[5,1], RC[5,2])$P
WK_compare_10 <- pressurecalculator_wk2(flow, time, RC[10,1], RC[10,2])$P
WK_compare_15 <- pressurecalculator_wk2(flow, time, RC[15,1], RC[15,2])$P
WK_compare_20 <- pressurecalculator_wk2(flow, time, RC[20,1], RC[20,2])$P

df_new_1 <- data.frame(R = rep(RC[1,1],length(flow_new)), C = rep(RC[1,2], length(flow_new)), Q = flow_new, t = time_new)
y_hat_1 <- predict(test_model, as.matrix(df_new_1),"SK")

df_new_5 <- data.frame(R = rep(RC[5,1],length(flow_new)), C = rep(RC[5,2], length(flow_new)), Q = flow_new, t = time_new)
y_hat_5 <- predict(test_model, as.matrix(df_new_5),"SK")

df_new_10 <- data.frame(R = rep(RC[10,1],length(flow_new)), C = rep(RC[10,2], length(flow_new)), Q = flow_new, t = time_new)
y_hat_10 <- predict(test_model, as.matrix(df_new_10),"SK")

df_new_15 <- data.frame(R = rep(RC[15,1],length(flow_new)), C = rep(RC[15,2], length(flow_new)), Q = flow_new, t = time_new)
y_hat_15 <- predict(test_model, as.matrix(df_new_15),"SK")

df_new_20 <- data.frame(R = rep(RC[20,1],length(flow_new)), C = rep(RC[20,2], length(flow_new)), Q = flow_new, t = time_new)
y_hat_20 <- predict(test_model, as.matrix(df_new_20),"SK")


#NB remember to change this x = time/time_sim when you change indices!!!!!!
df_plot <- data.frame(x = time_new, y1 = y_hat$mean,y2 = WK_compare, lb = y_hat$lower95, ub =y_hat$upper95, y3 = y_hat_1$mean, y4 = y_hat_5$mean, y5 = y_hat_10$mean, y9 = WK_compare_1, y10 = WK_compare_5, y11 = WK_compare_10, y12 = WK_compare_15, y13=WK_compare_20, y6 = y_hat_15$mean, y7 = y_hat_20$mean, y8 = (WK_compare_optimal-mu_WKcomp))

ggplot(df_plot,aes(x = Time(s),y = PredictedPressure))+
  geom_line(aes(x = x, y = y1, color = "My Surrogate Fit"))+
  geom_line(aes(x = x, y = y2, color = "WK2 model"))+
  #geom_line(aes(x = x, y = y8, color = "WK2 model-mean"))+
  geom_line(aes(x = x, y = y3), color = "grey", size = 0.3, linetype = "dashed")+
  geom_line(aes(x = x, y = y4), color = "grey", size = 0.3, linetype = "dashed")+
  geom_line(aes(x = x, y = y5), color = "grey", size = 0.3, linetype = "dashed")+
  geom_line(aes(x = x, y = y9), color = "#00CCCC")+
  geom_line(aes(x = x, y = y10), color = "#00CCCC")+
  geom_line(aes(x = x, y = y11), color = "#00CCCC")+
  
  geom_line(aes(x = x, y = y12), color = "#00CCCC")+
  geom_line(aes(x = x, y = y13), color = "#00CCCC")+
  geom_line(aes(x = x, y = y6), color = "grey", size = 0.3, linetype = "dashed")+
  geom_line(aes(x = x, y = y7), color = "grey", size = 0.3, linetype = "dashed")+
  
  geom_line(aes(x = x, y = ub), color = "#FF3366", linetype = "dashed")+
  geom_line(aes(x = x, y = lb), color = "#FF3366", linetype = "dashed")+
  theme(panel.background = element_rect(fill = NA, color = "black"))+
  ggsave("SurrogateFit.png")


cat("VARIANCE \n",test_model@covariance@sd2)
