#######################################################################################
#
#   Filename    :	      Table 2 and Table 5.R    												  
#                                                                                     
#   Project     :       Article "Estimating the correlation between semi-competing risk survival endpoints"                                                             
#   Authors     :       L Sorrell, Y Wei, M Wojtys and P Rowe                                                               
#   Date        :       01/06/2021
#																				  
#   R Version   :       R-3.6.1                                                              
#
#   Required R packages :  copula, mvtnorm, mstate
#
########################################################################################
library(copula) 
library(mvtnorm)
library(mstate) #ACS data
source("Functions.R")   #Likelihood functions
library(foreach)
data(aidssi2) #Load ACS data set

X <- aidssi2$si.time-aidssi2$entry.time #time to non-terminal event (SI switch)
d1 <- aidssi2$si.stat #indicator for SI switch
Y <- aidssi2$death.time-aidssi2$entry.time #time to terminal event (death from AIDS)
d2 <- aidssi2$death.stat #indicator for death from AIDS
df <- data.frame(X,Y,d1,d2)

########################
## Recreating Table 2 ## 
########################

table2 <- addmargins(table(d2,d1))

########################
## Recreating Table 5 ## 
########################

#prepare table as data frame:
lambda1 <- rep(NA, 4) #lambda_1
lambda1_lwci <- rep(NA, 4) #lower confidence interval for lambda_1
lambda1_upci <- rep(NA,4) #upper confidence interval for lambda_1
lambda2 <- rep(NA, 4) #lambda_2
lambda2_lwci <- rep(NA, 4) #lower confidence interval for lambda_2
lambda2_upci <- rep(NA,4) #upper confidence interval for lambda_2
rho <- rep(NA, 4) #rho
rho_lwci <- rep(NA, 4) #lower confidence interval for rho
rho_upci <- rep(NA,4) #upper confidence interval for rho
AIC <- rep(NA,4) #AIC
copula <- c("Normal", "Clayton", "Frank", "Gumbel")
get_initial_value<- function(para_list, likelihood_func, X, Y, d1, d2 ){

  combin_para= expand.grid(para_list)
  best_combin = 0
  temp_likelihood = -10000000
  for (variable in 1:nrow(combin_para)){
  # foreach(variable = 1:nrow(combin_para)) %dopar% {
    par = combin_para[variable,]

    likelihood = likelihood_func(par,X, Y, d1, d2)
      if ((!is.nan(likelihood))&&(is.finite(likelihood))){
         cat("likelihood:", likelihood, "\n")
         if(!is.na(temp_likelihood < likelihood)){
           if(temp_likelihood < likelihood){
             temp_likelihood <- likelihood
             best_combin = variable
             cat("best_combin:", best_combin, "\n") 
             cat("temp_likelihood:", temp_likelihood, "\n")
           }
         
        }
      }
  }
  cat("last_best_combin:", best_combin, "\n") 
  cat("last_temp_likelihood:", temp_likelihood, "\n")
  par_best = combin_para[best_combin,]
  return (par_best)
}

get_initial_value_parallel<- function(para_list, likelihood_func, X, Y, d1, d2 ){
  
  combin_para= expand.grid(para_list)
  # iterd <-iterators::iter(combin_para, by='row')
  a <-  foreach(para=iterators::iter(combin_para, by='row'), .combine='rbind') %dopar% {
    likelihood_func(para,X, Y, d1, d2)
  }
  value = max(a[is.finite(a)])
  
  a[!is.finite(a)] = -10000000000000
  value_index = which(a==max(a), arr.ind=TRUE)
  return (value_index)
}





get_initial_value_parallel1<- function(para_list, likelihood_func, X, Y, d1, d2 ){
  
  combin_para= expand.grid(para_list)
  # iterd <-iterators::iter(combin_para, by='row')
  a <-  foreach(para=iterators::iter(combin_para, by='row'), .combine='rbind') %dopar% {
    # likelihood_func(para,X, Y, d1, d2)
    l1 <- para[1]
    l2 <- para[2]
    theta <- para[3]
    
    C <- exp(-((-log(exp(-l1*X)))^(theta)+(-log(exp(-l2*Y)))^(theta))^(1/theta))
    C <- as.numeric(C)
    part1 <- ifelse(d1*d2==1,(log(C)+(theta-1)*log(-log(exp(-l1*X)))+(theta-1)*log(-log(exp(-l2*Y)))+log(theta-1+((-log(exp(-l1*X)))^theta+(-log(exp(-l2*Y)))^theta)^(1/theta))-log(exp(-l1*X))-log(exp(-l2*Y))-(2*theta-1)*log(-log(C))+log(l1)-l1*X+log(l2)-l2*Y),0)
    part2 <- ifelse(d1*(1-d2)==1,(log(C)+(theta-1)*log(-log(exp(-l1*X)))-log(exp(-l1*X))-(theta-1)*log(-log(C))+log(l1)-l1*X),0)
    part3 <- ifelse(((1-d1)*(d2))==1,(log(C)+(theta-1)*log(-log(exp(-l2*Y)))-log(exp(-l2*Y))-(theta-1)*log(-log(C))+log(l2)-l2*Y),0)
    part4 <- ifelse(((1-d1)*(1-d2))==1,log(C),0)
    
    loglik <- sum(as.numeric(part1)+as.numeric(part2)+as.numeric(part3)+as.numeric(part4) )#sum(part1+part2+part3+part4) 
  }
  value = max(a[is.finite(a)])
  
  a[!is.finite(a)] = -10000000000000
  value_index = which(a==max(a), arr.ind=TRUE)
  return (value_index)
}
par2<-par1 <- seq(-10 ,10,0.5)   #  0.01,10,1
par3<- seq(1,10,0.5) #. 0.1,1,0.1
para_list_gumbel= list(par1, par2, par3)
index_ini_gumbel = get_initial_value_parallel1(para_list_gumbel, gumbel_loglik, X, Y, d1, d2 ) 

# get_initial_value<- function(para_list, likelihood_func, X, Y, d1, d2 ){
# 
#   combin_para= expand.grid(para_list)
#   best_combin = 0
#   par_best=0
#   temp_likelihood = -10000000000
#   for (variable in 1:nrow(combin_para)) {
#     para = combin_para[variable,]
#     l1 <- para[1]
#     l2 <- para[2]
#     theta <- para[3]
#     
#     C <- -1/theta * log(((1-exp(-theta)-(1-exp(-theta*exp(-l1*X)))*(1-exp(-theta*exp(-l2*Y)))))/(1-exp(-theta)))
#     
#     part1 <- ifelse(d1*d2==1,(log(theta)+theta*C+log(exp(theta*C)-1)-log(exp(theta*exp(-l1*X))-1)-log(exp(theta*exp(-l2*Y))-1)+log(l1)-l1*X+log(l2)-l2*Y),0)
#     part2 <- ifelse(d1*(1-d2)==1,(log((1-exp(theta*C))/(1-exp(theta*exp(-l1*X))))+log(l1)-l1*X),0)
#     part3 <- ifelse(((1-d1)*(d2))==1,(log((1-exp(theta*C))/(1-exp(theta*exp(-l2*Y))))+log(l2)-l2*Y),0)
#     part4 <- ifelse(((1-d1)*(1-d2))==1,log(C),0)
#     
#     likelihood <- sum(as.numeric(part1)+as.numeric(part2)+as.numeric(part3)+as.numeric(part4) )#sum(part1+part2+part3+part4)
#     print(likelihood)
#     # likelihood = likelihood_func(par,X, Y, d1, d2)
#     if (!is.nan(likelihood))
#     {
#       is_samll = as.numeric(temp_likelihood) < as.numeric(likelihood)
#       #   cat("is_samll:", is_samll, "\n")
#       if(is_samll){
#         
#         temp_likelihood <- likelihood
#         best_combin <- variable
#         print(best_combin)
#         cat("temp_likelihood:", temp_likelihood, "\n")
#         
#       }
#     }
#   }
#   print(best_combin)
#   print(temp_likelihood)
# 
#   par_best = combin_para[best_combin,]
#   return (par_best)
# }

# par2<-par1 <- seq(0.1,5,0.2)   #  0.01,10,1
# par3<- seq(0.1,10,0.2) #. 0.1,1,0.1
# par2<-par1 <- seq(-10 ,10,0.5)   #  0.01,10,1
# par3<- seq(-10,10,0.5) #. 0.1,1,0.1
# para_list= list(par1, par2, par3)
# # para_list= list(par1, par2, par3, par1, par2)
# combin_para= expand.grid(para_list)
# 
# 
# best_par1 = get_initial_value_parallel(para_list, frank_loglik, X, Y, d1, d2 )  # clayton_loglik get_initial_value 
# 
# iterd <-iterators::iter(combin_para, by='row')
# a <-  foreach(para=iterators::iter(combin_para, by='row'), .combine='rbind') %dopar% {
#  frank_loglik(para,X, Y, d1, d2)
# }
# value = max(a[is.finite(a)])
# a[!is.finite(a)] = -10000000000000
# value_index = which(a==max(a),arr.ind=TRUE)




# 
# 
# 
# best_combin = 0
# 
# par = combin_para[23,]
# parxx= par[1]
# temp_likelihood = -100000
# likelihood = frank_loglik(par,X, Y, d1, d2)
# l1 <- par[1] #hazard rate for non-terminal event
# l2 <- par[2] #hazard rate for terminal event
# theta <- par[3] #association parameter
# 
# C <- (exp(-l1*X)^(-theta)+exp(-l2*Y)^(-theta)-1)^(-1/theta) #copula
# 
# part1 <- ifelse(d1*d2==1,(log(1+theta)+(1+2*theta)*log(C)-(theta+1)*log(exp(-l1*X))-(theta+1)*log(exp(-l2*Y))+log(l1)-l1*X+log(l2)-l2*Y),0) #both events
# part2 <- ifelse(d1*(1-d2)==1,((theta+1)*log(C)-(theta+1)*log(exp(-l1*X))+log(l1)-l1*X),0) #non-terminal event only
# part3 <- ifelse(((1-d1)*(d2))==1,((theta+1)*log(C)-(theta+1)*log(exp(-l2*Y))+log(l2)-l2*Y),0) #terminal event only
# part4 <- ifelse(((1-d1)*(1-d2))==1,log(C),0) #both events censored
# 
# likelihood <- sum(as.numeric(part1)+as.numeric(part2)+as.numeric(part3)+as.numeric(part4) )
# if(temp_likelihood < likelihood)
# {
#   
#   temp_likelihood = likelihood
#   best_combin = variable
#   print(best_combin)
#   print(temp_likelihood)
# }
#  
# 
# par_best = combin_para[best_combin,]

