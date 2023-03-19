# *************************************************** *********************************
#*********Qingmin Zhang writes this code ****************
#*********for utilising the coxPh model for data*********
#*
#*
#**********************************
# aidssi2, analysising the relationship between aid death and age and CCR5.

# install.packages(c("survival", "survminer", "gt", "conflicted","tidyverse","glue","webshot2"))

#library(ggplot2)
library(ggplot2)
library(ggpubr)
library("survminer")
library("survival")
library(mstate)
library(gt)


data("aidssi2")
print(aidssi2)
# To save the aidssi2 as .xls frormat
patnr = data.frame(aidssi2$patnr)
entry.time = data.frame(aidssi2$entry.time)
aids.time = data.frame(aidssi2$aids.time)
aids.stat = data.frame(aidssi2$aids.stat)
si.time = data.frame(aidssi2$si.time)
si.stat = data.frame(aidssi2$si.stat)
death.time = data.frame(aidssi2$death.time)
death.stat = data.frame(aidssi2$death.stat)
age.inf = data.frame(aidssi2$age.inf)
ccr5 = data.frame(aidssi2$ccr5)
aidssi2 = cbind(patnr,
                entry.time,
                aids.time,
                aids.stat,
                si.time,
                si.stat,
                death.time,
                death.stat,
                age.inf,
                ccr5)
names(aidssi2) <- c("patnr",
                    "entry.time",
                    "aids.time",
                    "aids.stat",
                    "si.time",
                    "si.stat",
                    "death.time",
                    "death.stat",
                    "age.inf",
                    "ccr5")
write.csv(data.frame(aidssi2), file ="aidssi2.csv")
# discri_statistic = data.frame(summary(aidssi2))

# kaplan-meier curve for original data
Y <- data.frame(aidssi2$death.time-aidssi2$entry.time) #time to terminal event (death from AIDS)
d2 <- data.frame(aidssi2$death.stat)#indicator for death from AIDS
age<- data.frame(aidssi2$age.inf)
ccr5 <- data.frame(aidssi2$ccr5 )
aidssi2_frame = cbind(Y, d2, age, ccr5)
write.csv(data.frame(aidssi2_frame), file ="aidssi2_frame.csv")
names(aidssi2_frame) <- c("Y", "d2", "age","ccr5")
km<-survfit(Surv(Y, d2)~1, data=aidssi2_frame)
# set.seed(33)
# palette <- sample(c("color1", "color2", ...), 324, replace = TRUE)
ggsurvplot(km, 
           title  = "kaplan-meier curve",  
           xlab = "time(year)", 
           ylab = "survival probability", 
           font.main = c(16, "bold", "darkblue"), 
           font.x = c(12,  "bold.italic","darkblue"), 
           font.y = c(14, "darkblue"),
           legend.labs = ("all"),  
           size = 1, 
           break.x.by=12 ,
           break.y.by=0.2 ,
           surv.scale="percent" ,
           # palette ="aaas", 
           conf.int = T, 
           pval = FALSE ,
           pval.coord = c(12, 0),
           pval.size = 5,
           pval.method=TRUE,
           pval.method.size=5,
           pval.method.coord=c(1,0),
           surv.median.line = "hv",
           ggtheme = theme_bw(), 
           palette= "aaas", 
           ylim = c(0, 1)
)
####### plot the predicted survival function curve
# explicitly change the dummy variables
Y <- data.frame(aidssi2$death.time-aidssi2$entry.time) #time to terminal event (death from AIDS)
d2 <- data.frame(aidssi2$death.stat)#indicator for death from AIDS
age<- data.frame(aidssi2$age.inf)
ccr5 <- data.frame(aidssi2$ccr5 )
aidssi2_frame = cbind(Y, d2, age, ccr5)
names(aidssi2_frame) <- c("Y", "d2", "age","ccr5")
# kaplan-meier curve for original data
km<-survfit(Surv(Y, d2)~1, data=aidssi2_frame)
# set.seed(33)
# palette <- sample(c("color1", "color2", ...), 324, replace = TRUE)
ggsurvplot(km, 
           #title  = "kaplan-meier curve",  
           ggtheme = theme_minimal()
)
####### plot the predicted survival function curve
# explicitly change the dummy variables
modelData <- data.frame(model.matrix(~ccr5+Y+d2+ age  , aidssi2_frame)  )
names(modelData) <- c("intercept","ccr5", "Y", "d2", "age")
# model.matrix(~ccr5, aidssi2)
# coxph model 
res.cox <- coxph(Surv(Y, d2) ~ age + ccr5, data =  modelData)
res_c0x_sum_dummy = summary(res.cox)
ggsurvplot(survfit(res.cox), data =  modelData, palette= "#2E9FDF",
           ggtheme = theme_minimal())
# implicitly deal with the dummy variables 
res.cox <- coxph(Surv(Y, d2) ~ age + ccr5, data =  aidssi2_frame)
res_c0x_sum = summary(res.cox)

ggsurvplot(survfit(res.cox), data =  aidssi2_frame,palette= "#2E9FDF",
           ggtheme = theme_minimal())

# %>%
# dplyr::slice(1:5)
sex_df <- with(aidssi2_frame,data.frame(ccr5 = c('WW', 'WM'), 
                                        age = rep(mean(age, na.rm = TRUE), 2)))
sex_df

fit_rescox <- survfit(res.cox, newdata = sex_df)
ggsurvplot(fit_rescox, data =sex_df, conf.int = TRUE, legend.labs=c("ccr5=WW", "ccr5=WM"),
           ggtheme = theme_minimal())
# save the estimation of the coefficient
dataframe_res_c0x = data.frame(res_c0x_sum[["coefficients"]])
dataframe_res_cox_sum = data.frame(res_c0x_sum[["conf.int"]])
bind_cox = cbind(round(dataframe_res_c0x,3),  round(dataframe_res_cox_sum,3) )
bind_cox1= bind_cox[,-which(duplicated(names(bind_cox)))]
write.csv(data.frame(bind_cox1), file ="coefficient_Cox.csv")
#
logtest =  res_c0x_sum[["logtest"]] 
sctest = res_c0x_sum["sctest"]
waldtest = res_c0x_sum["waldtest"]
test_combin = data.frame( logtest,  sctest , waldtest  )
test_combin = data.frame(t(test_combin))

