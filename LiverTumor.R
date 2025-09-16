library(readxl)
library(tidyverse)
library(modelsummary)
library(lattice)
library(survival)
library(ggfortify)
library(eventglm)
library(timereg)
library(cmprsk)
library(Hmisc)

liver <- read_excel("D:/Uni/KU - Courses/Advanced Topics in Survival Analysis/Liver Tumor.xlsx")
relapse <- read_excel("D:/Uni/KU - Courses/Advanced Topics in Survival Analysis/Liver Tumor.xlsx",sheet = "Foglio2")

# Data Exploration ####
glimpse(liver,width = 62)
table(liver$deathByTumor)

liver$tumorDimension <- as.numeric(liver$tumorDimension)




status <- liver$deathByTumor+liver$death
# 0 = censored
# 1 = Dead
# 2 = Death by Tumor
table(status)
group <- factor(liver$tumorBack) 


## Conditional Distribution ####
liver %>% 
  mutate(status = status) %>% 
  group_by(tumorBack, status) %>% 
  summarize_all(mean, na.rm = T) %>% 
  select(status, tumorDimension)

liver %>% 
  mutate(
    fstatus = factor(case_when(
      status == 0 ~ "Censored",
      status == 1 ~ "DOC",
      TRUE~ "Tumor Death"
    )),
    postSurgeryComplications = factor(ifelse(postSurgeryComplications==1, "Had Post Surgery Complications", "No Complications")),
    bloodLossDuringSurgery = factor(ifelse(bloodLossDuringSurgery==1, "Had Blood Loss During Surgery", "No Blood Loss"))
    ) %>% 
  ggplot(aes(x = lengthSurgery, group = fstatus))+
  geom_boxplot(stat = "boxplot", varwidth =T , aes(fill = fstatus), alpha = 0.2)+
  coord_flip()+
  scale_y_discrete() +
  facet_wrap(~postSurgeryComplications + bloodLossDuringSurgery,ncol = 4, nrow = 1)+
  theme_bw()

liver %>% 
  mutate(
    fstatus = factor(case_when(
      status == 0 ~ "Censored",
      status == 1 ~ "DOC",
      TRUE~ "Tumor Death"
    )),
    hypertension = factor(ifelse(hypertension==1, "Has Hypertension", "Does Not Have Hypertension")),
    ) %>% 
  ggplot(aes(x = age, group = fstatus))+
  geom_boxplot(stat = "boxplot", varwidth = F, aes(fill = fstatus), alpha = 0.2)+
  coord_flip()+
  scale_y_discrete() +
  facet_wrap(~hypertension,ncol = 2, nrow = 1)+
  theme_bw()

liver %>% 
  mutate(
    fstatus = factor(case_when(
      status == 0 ~ "Censored",
      status == 1 ~ "DOC",
      TRUE~ "Tumor Death"
    )),
    tumorBack = factor(ifelse(tumorBack==1, "Tumor Relapsed", "Tumor Did Not Relapsed")),
    ) %>% 
  ggplot(aes(x = tumorDimension, group = fstatus))+
  geom_boxplot(stat = "boxplot", varwidth = T, aes(fill = fstatus), alpha = 0.2)+
  coord_flip()+
  scale_y_discrete() +
  xlim(c(0,100))+
  facet_wrap(~tumorBack,ncol = 2, nrow = 1)+
  theme_bw()


bwplot(age~status | hypertension , data = liver, horizontal = F)
bwplot(lengthSurgery~status | postSurgeryComplications , data = liver, horizontal = F)
bwplot(tumorDimension~status | tumorBack , data = liver, horizontal = F)
par(mfrow = c(1,1))

## Pairs Plot ####

panel.smooth <- function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
          cex = 1, col.smooth = 1:3, span = 2/3, iter = 3, ...){
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) 
    lines(stats::lowess(x[ok & status ==0], y[ok & status == 0], f = span, iter = iter), 
          col = col.smooth[1], ...)
    lines(stats::lowess(x[ok & status ==1], y[ok & status == 1], f = span, iter = iter), 
          col = col.smooth[2], ...)
    lines(stats::lowess(x[ok & status ==2], y[ok & status == 2], f = span, iter = iter), 
          col = col.smooth[3], ...)
}


panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y,use = "complete.obs"))
  sign <- ifelse(cor(x, y,use = "complete.obs")>0,"","-")
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(sign,prefix, txt)
  cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

panel.dens <- function(x,...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5))
  
  # Compute the density
  dens <- density(na.omit(x))
  
  # Normalize the density for plotting
  y <- dens$y / max(dens$y)
  
  # Plot the density
  lines(dens$x, y, color = 1)
}



pairs(liver[,c("daysSinceSurgery","age","tumorDimension", "lengthSurgery")],
      col = status+1,
      lwd = 2,
      pch = "*",
      lower.panel = panel.smooth,
      upper.panel= panel.cor,
      diag.panel = panel.dens
      )
legend("center",pch = "*" ,col  = 1:3,legend = c("Censored", "DOC", "Tumor Death"), bty = "n")


## Cure model ####

cure <- survfit(Surv(relapse$daysToRelapse,relapse$Relapse)~1)
plot(cure,mark.time=T,conf.int=F)

# Cumulative Incidence Function #####

status1 <- ifelse(status != 1, 0,1)
status2 <- ifelse(status != 2, 0,1)

KM0 <- survfit(Surv(daysSinceSurgery,status1)~1,data=liver, type="kaplan-meier")
KM1 <- survfit(Surv(daysSinceSurgery,status2)~1,data=liver, type="kaplan-meier")

fitMARG <- cuminc(ftime = liver$daysSinceSurgery,
              fstatus = status,
              cencode = 0)

par(mfrow = c(1,2))
plot(fitMARG$`1 2`$time, fitMARG$`1 2`$est,col=2,type="l",main="Death by Tumor",xlab="Time", ylab="Probability", ylim =c(0,1))
lines(KM1$time, 1- KM1$surv)
legend("topleft",legend=c("Marginal Distribution","Sub Distribution"),lty = 1, col=1:2, bty = "n", cex = .8)
plot(fitMARG$`1 1`$time, fitMARG$`1 1`$est,col=2,type="l",main="Other Causes",xlab="Time", ylab="Probability", ylim =c(0,1))
lines(KM0$time, 1- KM0$surv)
par(mfrow = c(1,1))





fit <- cuminc(ftime = liver$daysSinceSurgery,
              fstatus = status,
              group = liver$tumorBack,
              cencode = 0)

est.NOTB.KOTHER <-  fit$`0 1`$est
time0 <- fit$`0 1`$time

est.TB.KOTHER <- fit$`1 1`$est
time1 <- fit$`1 1`$time

est.NOTB.KT <- fit$`0 2`$est
time2 <- fit$`0 2`$time

est.TB.KT <- fit$`1 2`$est
time3 <- fit$`1 2`$time


plot(time2,est.NOTB.KT,col=2, lwd = 2,type="l",main="Cumulative Incidence Function",xlab="Time", ylab="Probability", ylim =c(0,.6))
lines(time0,est.NOTB.KOTHER, lwd = 2,col=1)
lines(time1,est.TB.KOTHER,col=1, lwd = 2,lty = "dashed")
lines(time3,est.TB.KT,col=2, lwd = 2,lty = "dashed")
legend("topleft",lwd = 2,legend=c("Tumor Death No Relapse","Other Cause No Replace", "Tumor Death with Relapse", "Other Cause with Relapse"), col=c(2,1,2,1),lty=c("solid","solid","dashed", "dashed"), bty = "n", cex = 1)

liver$id <- 1:nrow(liver)
liver$status <- status
colSums(is.na(liver))
# na.liver <- na.omit(liver[,c(1,4:7,9:14)])
na.liver <- na.omit(liver[,-c(2,3,4,8,9)])
mod <- na.liver$lengthSurgery%/%60
na.liver$lengthSurgery <- mod+(na.liver$lengthSurgery-mod*60)/60
# na.liver <- na.omit(liver[,-c(2,3)])

# Regression with simple COX PH ####

mdl.ph1 <- coxph(Surv(na.liver$daysSinceSurgery,status[na.liver$id]==1)~.,control =  list(iter.max = 1e4,timefix=T), data = na.liver[,-9])
mdl.ph2 <- coxph(Surv(na.liver$daysSinceSurgery,status[na.liver$id]==2)~.,control =  list(iter.max = 1e4,timefix=T), data = na.liver[,-9])


(sph1 <- summary(mdl.ph1))
(sph2 <- summary(mdl.ph2))
## stepAIC() ####

step1 <- MASS::stepAIC(mdl.ph1, direction = "both",)
step2 <- MASS::stepAIC(mdl.ph2, direction = "both")
# step1$anova # lengthSurgery
# step2$anova # sex + hypertension + lengthSurgery + tumorBack
# 
# 
# 
# ## Results ####
# ph1 <- coxph(Surv(na.liver$daysSinceSurgery,status[na.liver$id]==1)~ lengthSurgery,control =  list(iter.max = 1e4,timefix=T), data = na.liver, model = T)
# ph2 <- coxph(Surv(na.liver$daysSinceSurgery,status[na.liver$id]==2)~ sex + hypertension + lengthSurgery + tumorBack,control =  list(iter.max = 1e4,timefix=T), data = na.liver)
# 
# summary(ph1)
# summary(ph2)
# 
# plot(ph1$y, add = T, main = "Predicted Survival Function")
# par(new=TRUE)
# plot(ph2$y, col = 2)
# legend("bottomleft", col = 1:2, lty = c(1,1), legend = c("DOC", "Tumor Death"), bty = "n" )
# 

# Regression on Sub Distribution Hazards ####
X <- na.liver[,2:8]
Y <- na.liver$daysSinceSurgery
D <- status[na.liver$id]

fit1 <- crr(Y,D,X,failcode=1,cencode = 0, maxiter = 1e3)
fit2 <- crr(Y,D,X,failcode=2,cencode = 0, maxiter = 1e3)
(s1 <- summary(fit1))
(s2 <- summary(fit2))



## Survival curve for a certain profile #####
X1=c(65,0,1,7,1,1,1)
X2=c(65,0,1,7,1,1,0)
X3=c(65,0,1,2,1,1,0)
X4=c(65,0,1,2,1,1,1)


pc1 <- predict(fit1,rbind(X1,X2,X3,X4))
pc2 <- predict(fit2,rbind(X1,X2,X3,X4))

par(mfrow=c(1,2))
plot(pc1,col=c(1:4),lty = 1,lwd = 2, main="Non-Cancer Death",xlab="Time",ylab="Probability", ylim = c(0,1))
legend("topleft",bty = "n",legend=c("7 Hours Surgery + TumorBack", "7 Hours Surgery + No TumorBack","2 Hours Surgery + No TumorBack", "2 Hours Surgery + TumorBack"),cex = .8, col=1:4,lty=1)
plot(pc2,col=1:4,lty = 1,lwd = 2, main="Cancer Death",xlab="Time", ylab="Probability", ylim = c(0,1))


# Proportional Odds Model ####
odds.subd1 <- prop.odds.subdist(Event(daysSinceSurgery, status)~age + sex+ hypertension + lengthSurgery + bloodLossDuringSurgery+postSurgeryComplications+ tumorBack,
                  data = na.liver, cause = 1)

odds.subd2 <- prop.odds.subdist(Event(daysSinceSurgery, status)~age + sex+ hypertension + lengthSurgery + bloodLossDuringSurgery+postSurgeryComplications+ tumorBack,
                  data = na.liver, cause = 2)

(sd1 <- odds.subd1)
(sd2 <- odds.subd2)

newX <- rbind(X1,X2,X3,X4)
colnames(newX) <- names(fit1$coef)
psd1 <- predict(odds.subd1, Z = newX)
psd2 <- predict(odds.subd2, Z = newX)

par(mfrow=c(1,2))
matplot(x = psd1$time, t(psd1$P1), col = 1:4, type = "l", lty = 1, main="Non-Cancer Death",xlab="Time", ylab="Probability", ylim = c(0,1))
legend("topleft",bty = "n",legend=c("7 Hours Surgery + TumorBack", 
                                    "7 Hours Surgery + No TumorBack","2 Hours Surgery + No TumorBack", "2 Hours Surgery + TumorBack"),
       cex = .8, col=1:4,lty=1)
matplot(x = psd2$time, t(psd2$P1), col = 1:4, type = "l", lty = 1, main="Cancer Death",xlab="Time", ylab="Probability", ylim = c(0,1))


# Total Results ####

xx <- as.data.frame(
cbind(
  covariate = rownames(sph1$coef),
  coeff1ph = paste0(round(sph1$coef[,2],3), gtools::stars.pval(sph1$coef[,5])),#, " (", round(sph1$coef[,3],3),")"),
  coeff2ph = paste0(round(sph2$coef[,2],3), gtools::stars.pval(sph2$coef[,5])),#, " (", round(sph2$coef[,3],3),")"),
  coeff1sd = paste0(round(s1$coef[,2],3), gtools::stars.pval(s1$coef[,5])),#,"(", round(s1$coef[,3],3),")"),
  coeff2sd = paste0(round(s2$coef[,2],3), gtools::stars.pval(s2$coef[,5])),#,  " (", round(s2$coef[,3],3),")"),
  coeff1po = paste0(round(exp(coef(sd1)[,1]),3), gtools::stars.pval(coef(sd1)[,6])),#, " (", round(coef(sd1)[,2],3),")"),
  coeff2po = paste0(round(exp(coef(sd2)[,1]),3), gtools::stars.pval(coef(sd2)[,6]))#, " (", round(coef(sd2)[,2],3),")")
      )
)

rownames(xx) <- xx$covariate
xx <- xx[,-1]

xtable::xtable((xx))

