library(rms) #RCS
library(survminer)#曲线
library(ggplot2)#画图
library(ggsci) #调色板
library(mstate)
library(gt)
dt <- data("aidssi2")
patnr = data.frame(aidssi2$patnr)
entry.time = data.frame(aidssi2$entry.time)
aids.time = data.frame(aidssi2$aids.time)
aids.stat = data.frame(aidssi2$aids.stat)

si.time = data.frame(aidssi2$si.time)
si.stat = data.frame(aidssi2$si.stat)
death.time = data.frame(aidssi2$death.time)
death.stat = data.frame(aidssi2$death.stat)
age.inf = data.frame(aidssi2$age.inf)
# set the analysis data
Y <- data.frame(aidssi2$death.time+aidssi2$entry.time) #infection time to terminal event (death from AIDS)
d2 <- data.frame(aidssi2$death.stat)#indicator for death from AIDS
age<- data.frame(aidssi2$age.inf)
ccr5 <- data.frame(aidssi2$ccr5 )
aidssi2_frame = cbind(Y, d2, age, ccr5)
write.csv(data.frame(aidssi2_frame), file ="aidssi2_frame.csv")
names(aidssi2_frame) <- c("time", "state", "age","ccr5")


x_start = 0
x_end = max(aidssi2_frame$time)+0.5

#设定数据环境
dd<-datadist(aidssi2_frame) 
options(datadist='dd')

S <- Surv(aidssi2_frame$time, aidssi2_frame$state)
Hb = aidssi2_frame$age
ccr5 = data.frame(aidssi2$ccr5)
##构建模型
fit <- cph(S ~ rcs(aidssi2_frame$age,4) + ccr5, x=TRUE, y=TRUE,data=aidssi2_frame)
cox.zph(fit, "rank")             # PH 检验
ggcoxzph(cox.zph(fit, "rank"))   #可视化等比例假定
anova(fit)                       #非线性检验



Pre0 <-rms::Predict(fit,rcs(aidssi2_frame$age,4) ,fun=exp,type="predictions",ref.zero=TRUE,conf.int = 0.95,digits=2);
##其中fun是转化函数
ggplot(Pre0)
View(Pre0)

ggplot()+
  geom_line(data=Pre0,
            aes(Hb,yhat,colour=gender),alpha=0.7)+
  scale_color_nejm()+ ##采用ggsci包中英格兰调色，也可以其他
  geom_ribbon(data=Pre0,
              aes(Hb, ymin=lower,ymax=upper,fill=gender),alpha=0.1)+
  scale_fill_nejm()+
  geom_hline(yintercept=1,linetype=2,size=0.75)


Pre1 <- rms::Predict(fit, Hb, gender=c('male','female'),fun=exp,type="predictions",ref.zero=TRUE,conf.int = 0.95,digits=2)
par(mfrow=c(1,2))
ggplot(Pre1)
View(Pre1)