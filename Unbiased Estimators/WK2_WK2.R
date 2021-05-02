library(ggplot2)
library(reticulate)
np = import("numpy")

#setwd("Z:/PROSJEKTOPPGAVE/Script")

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

create_field_wk2 = function(flow, time, Rreal, Creal, sdnoise = 3, cycles, ind, seed = 123){
  Pwk2 = pressurecalculator_wk2(flow,time,Rreal,Creal)$P
  l = cycles*length(ind)
  df = data.frame(y = rep(NA, l), t = rep(time[ind],cycles), Q = rep(flow[ind], cycles))
  set.seed(seed)
  y = matrix(NA, nrow = length(ind), ncol = cycles)
  for(i in 1:cycles)  y[,i] = Pwk2[ind] + rnorm(length(ind), 0, sdnoise)
  df$y = as.vector(y)
  return(list(field_data = df, Pwk2 = data.frame(t = time,P = Pwk2)))
}

Sims2to2 <- function(flow, time,R,C, cycles2, sdnoise2, ind2, model){
  min.RSS <- function(data, par2){
    with(data, sum((p-pressurecalculator_wk2(flow,time,par2[1],par2[2])$P[ind2])^2))
  }
  #Store real values for plot later
  R_real <- R
  C_real <- C
  #initial values
  par2 <- c(0.8,2) #Initial values for R,C
  
  cycles2 <-cycles2 
  sdnoise2 <- sdnoise2 
  ind2 <- ind2 #indices of data points chosen to be same length as data to avoid non-conformable arrays in pressureCalc
  
  field_mat2 <- matrix(NA, nrow = cycles2, ncol = length(ind2)) #matrix of nr sims x observations 
  P_class2 <- rep(NA, cycles2) #pressure vector
  
  WK_real <- create_field_wk2(flow,time,R,C, sdnoise2, cycles2, ind2)
  np$max(WK_real$field_data$y)
  np$min(WK_real$field_data$y)
  
  WK2_run <- WK_real
  
  #sort field observations in ind times cycles matrix
  y_mat2 <- matrix(WK2_run$field_data$y,nrow = length(ind2), ncol = cycles2)
  
  #Find parameters for each cycle
  out_mat2 <- matrix(NA,nrow =cycles2, ncol = length(par2))
  
  #Sims
  WK_sims <- matrix(NA, nrow=cycles2,ncol=101)
  b_hat <- matrix(NA, nrow=cycles2,ncol = length(ind2))
  
  
  for (k in 1:cycles2){
    dat_run2 <- data.frame(p = y_mat2[,k])#sample diff seeds field obs
    field_mat2[k,] <- y_mat2[,k] #add field obs i to matrix for later storage
    
    
    P_class2[k] <- np$max(y_mat2[,k]) - np$min(y_mat2[,k])
    cat("Max P: \ ", np$max(y_mat2[,k])," \ ")
    cat("Min P: \ ", np$min(y_mat2[,k]), "cycle \ ", k,"\n" )
    
    res2 <- optim(par2, min.RSS, data = dat_run2) #optimize using initial parameters, and current field obs 
    
    for (j in 1:length(par2))out_mat2[k,j] = res2$par[j]
    
    #Save Pressure estimate for later
    P_cycle <- pressurecalculator_wk2(flow,time,res2$par[1],res2$par[2])$P
    WK_sims[k,] <- P_cycle
    
    #Calculate bias
    b_hat[k,] <- dat_run2$p-P_cycle[ind2]
    
    
  }
  df_R2 <- data.frame(R = out_mat2[,1])
  df_C2 <- data.frame(C = out_mat2[,2])
  
  R_mean2 <- sum(out_mat2[,1])/cycles2
  C_mean2 <- sum(out_mat2[,2])/cycles2
  
  ggplot(data=df_R2, aes(R)) + 
    geom_histogram()+
    geom_vline(aes(xintercept = R_real, color = "True R"))+
    geom_vline(aes(xintercept = R_mean2, color = "Estimated R"))+
    theme(panel.background = element_rect(fill = NA, color = "black"))+
    ggsave(stringr::str_c(model,"histogramR2to2.png"))
  
  ggplot(data=df_C2, aes(C)) + 
    geom_histogram()+
    geom_vline(aes(xintercept = C_real, color = "True C"))+
    geom_vline(aes(xintercept = C_mean2, color = "Estimated C"))+
    theme(panel.background = element_rect(fill = NA, color = "black"))+
    ggsave(stringr::str_c(model,"histogramC2to2.png"))
  
  #Boxplots 
  par_mat <- data.frame(R = out_mat2[,1],C = out_mat2[,2])
  
  ggplot(data = par_mat, aes(y = value, color = parameter))+
    geom_boxplot(aes(y = R, color = "R"))+
    theme(panel.background = element_rect(fill = NA, color = "black"))+
    ggsave(stringr::str_c(model,"boxR2to2.png"))
  
  ggplot(data = par_mat, aes(y = value, color = parameter))+
    geom_boxplot(aes(y = C, color = "C"))+
    theme(panel.background = element_rect(fill = NA, color = "black"))+
    ggsave(stringr::str_c(model,"boxC2to2.png"))
  
  #PI
  PI_WK_sims <- apply(WK_sims,2,sort)
  lower <- PI_WK_sims[5,]
  upper <- PI_WK_sims[95,]
  
  #Calculating final model profiles
  optimal_P <- pressurecalculator_wk2(flow,time,R_mean2,C_mean2)
  true_P <- pressurecalculator_wk2(flow,time,R,C)
  
  #mean bias 
  b_mean <- apply(abs(b_hat),2,sum)/100
  #plotting bias 
  df_bias <- data.frame(x = time[ind2],y1 = b_hat[100,],y2 = b_mean, y3 = b_hat[1,], y4 = b_hat[25,],y5 = b_hat[50,], y6 = b_hat[75,])
  
  
  df_field <- data.frame(x = WK_real$field_data$t, y1 = WK_real$field_data$y)
  df_optim <- data.frame(x = time, y1 = optimal_P$P, y2  = true_P$P, lb = lower, ub = upper)
  
  ggplot(df_field, aes(x = Time(s), y = AorticPressure(mmHg)))+
    geom_point(aes(x = x, y =y1),color = "grey")+
    geom_line(aes(x=x, y=y1,color = "Predicted Mean P"), data = df_optim)+
    geom_line(aes(x=x, y=y2,color = "Real P"), data = df_optim)+
    geom_line(data = df_optim, aes(x = x, y = lb), col ="#33CC66", size = 1, linetype = "dashed") + 
    geom_line(data = df_optim, aes(x = x, y = ub), col ="#33CC66", size = 1, linetype = "dashed")+
    geom_line(data = df_bias, aes(x = x, y = y1),color = "grey",  size = 1, linetype = "dashed")+
    geom_line(data = df_bias, aes(x = x, y = y3), color = "grey", size = 1, linetype = "dashed")+
    geom_line(data = df_bias, aes(x = x, y = y4), color = "grey", size = 1, linetype = "dashed")+
    geom_line(data = df_bias, aes(x = x, y = y5), color = "grey", size = 1, linetype = "dashed")+
    geom_line(data = df_bias, aes(x = x, y = y6), color = "grey", size = 1, linetype = "dashed")+
    geom_line(data = df_bias, aes(x = x, y = y2, col = "Mean Bias"), size = 1, linetype = "dashed")+
    theme(panel.background = element_rect(fill = NA, color = "black"))+
    ggsave(stringr::str_c(model,"fit2to2.png"))
  
  return(list(out_mat2 = out_mat2 ,R_mean2 = R_mean2, C_mean2 = C_mean2))
  
  
}

Flow_dat <- readRDS("WK2DatNew.rds")
flow <- Flow_dat[,2]*0.7
time <- Flow_dat[,1]

plot(time,flow,xlab = "time(s)",ylab="flow(ml/min)")

#Run for HP and nonHP
R <- 1.4
C_HP <- 0.8
C_nonHP <- 1.7

cycles2 <- 100
sdnoise2 <- 4
ind2 <- round(seq(1,101,length.out = 20))

WK_test <- create_field_wk2(flow,time,R,C_HP, sdnoise2, cycles2, ind2)
np$max(WK_test$field_data$y)
np$min(WK_test$field_data$y)

HP <- Sims2to2(flow,time,R,C_HP,cycles2,sdnoise2,ind2,"HP")
non_HP <- Sims2to2(flow,time,R,C_nonHP,cycles2,sdnoise2,ind2, "non_HP")


#Compare Compliance for HP and no HP 
par_matC_HP <- data.frame(variable = rep("C - HT", cycles2), value=HP$out_mat2[,2])
par_matC_HP
par_matC_nonHP <- data.frame(variable = rep("C - nonHT", cycles2), value= non_HP$out_mat[,2])
par_matC_nonHP

par_ultimateC <- rbind(par_matC_HP,par_matC_nonHP)
par_ultimateC

ggplot(par_ultimateC, aes(x=variable, y = value))+
  geom_boxplot()+
  ggsave("Cdiff2to2.png")



