
#Data
ying1 <- read.csv('https://raw.githubusercontent.com/VNyaga/ESMARConf2023/main/Ying1.csv', sep=',')

#metafor
#============================================================================================
library(metafor)
ying1$dis = with(ying1, tp + fn)
ying1$ndis = with(ying1, tn + fp)

se = rma.glmm(xi=tp, ni=dis, data=ying1, measure ="PLO")
sp = rma.glmm(xi=tn, ni=ndis, data=ying1, measure ="PLO")
sp
forest.rma(sp)

#0-1 scale
c(exp(sp$b)/(1 + exp(sp$b)), exp(sp$ci.lb)/(1 + exp(sp$ci.lb)), exp(sp$ci.ub)/(1 + exp(sp$ci.ub)))


se = rma.glmm(xi=tp, ni=dis, data=ying1, measure ="PLO", method="EE")
se

forest.rma(se)
#0-1 scale
c(exp(se$b)/(1 + exp(se$b)), exp(se$ci.lb)/(1 + exp(se$ci.lb)), exp(se$ci.ub)/(1 + exp(se$ci.ub)))


#MADA
#============================================================================================

library(mada)

library(stringr)
names(ying1)[2:5] <- str_to_upper(names(ying1)[2:5])
fit = reitsma(ying1, correction = 0)


windows()
par(mfcol=c(1,2))
forest(madad(ying1), type = "sens", snames=ying1$author, xlab="Sens")
forest(madad(ying1), type = "spec", snames=ying1$author, xlab="Spec")

fit = reitsma(ying1)
#0-1 scale
c(ss$coefficients[3,1], ss$coefficients[3,5], ss$coefficients[3,6])
c(1-ss$coefficients[4,1], 1-ss$coefficients[4,6], 1-ss$coefficients[4,5])


ss = summary(fit)
windows()
plot(fit, sroclwd = 2)
points(fpr(ying1), sens(ying1), pch = 2)
legend("bottomright", c("data", "summary estimate"), pch = c(2,1))
legend("bottomleft", c("SROC", "conf. region"), lwd = c(2,1))




#============================================================================================
df = read.csv('C:/DATA/WIV/Projects/Stata/Metadta/Data/dementia.csv')


fit = reitsma(df)

summary(fit)
windows()
plot(fit, sroclwd = 2)
points(fpr(df), sens(df), pch = 2)
legend("bottomright", c("data", "summary estimate"), pch = c(2,1))
legend("bottomleft", c("SROC", "conf. region"), lwd = c(2,1))

#============================================================================================
#Meta-regression
df = read.csv('C:/DATA/WIV/Projects/Stata/Metadta/Data/ascus.csv')
(fit <- reitsma(df, formula = cbind(tsens, tfpr) ~ Test))
summary(fit)


#============================================================================================
#Meta-regression
library(stringr)
df = read.csv('C:/DATA/WIV/Projects/Stata/Metadta/Data/midas_example.csv')
names(df)[3:6] <- str_to_upper(names(df)[3:6])

#Only 1 covariate allowed in mada

(fit <- reitsma(df, formula = cbind(tsens, tfpr) ~ prodesign + age))

(fit <- reitsma(df, formula = cbind(tsens, tfpr) ~ prodesign + age + qscore + sampsize + fulverif + testdescr + refdescr + subjdescr + report + brdspect + blinded))
summary(fit)


plot(sroc(fit, type = "ruttergatsonis"), ylim = c(0, 1) )

plot(x = df$fpr[df$Test=="HC2"], y = df$se[df$Test=="HC2"], col='red', xlim = c(0, 1), ylim = c(0, 1))
points(x = df$fpr[df$Test=="RepC"], y = df$se[df$Test=="RepC"], col='blue')

points(y=mu[1], x = mu[2], pch=19, col='red')
points(y=mu[1], x = mu[2], pch=19, col='blue')

lines(x=Uniroc$fpr, y =Uniroc$se, col="green")

#
ROCellipse(fit, add=TRUE)

lines(x = fpr, y = se, col='green')

lines(x = pfpr, y = pse, col='green')

#Telomerase
plot(x = df$fpr, y = df$se, col='black', xlim = c(0, 1), ylim = c(0, 1))

plot(x = df$se, y = df$sp, col='black', xlim = c(0, 1), ylim = c(0, 1))

plot(y = sesp$se, x = sesp$sp, col='purple', xlim = c(0, 1), ylim = c(0, 1), pch=19,
     xlab='Specificity', ylab = 'Sensitivity', cex.axis=1.5, cex.lab=1.5)

#plot(y = sesp$logitse, x = sesp$logitsp, col='black',  pch=19,xlim = c(-1, 5), ylim = c(-1, 5),
#     xlab='Logit Specificity', ylab = 'Logit Sensitivity', cex.axis=1.5, cex.lab=1.5)

points(y=0.75, x = 0.86, pch=19, col='red')

points(y=0.76, x = 0.89, pch=19, col='blue')

points(y=0.76, x = 0.81, pch=19, col='green')

arrows(0.8, 0.9, 0.86, 0.75, length=0.1)
text(0.78, 0.92, labels = 'Glas')

arrows(0.8, 0.5, 0.89, 0.76, length=0.1)
text(0.8, 0.48, labels = 'Riley')

arrows(0.6, 0.7, 0.81, 0.76, length=0.1)
text(0.55, 0.7, labels = 'Nyaga')

arrows(0.8, 0.5, 0.89, 0.76, length=0.1)
text(0.8, 0.48, labels = 'Riley')
lines(x=Uniroc$fpr, y =Uniroc$se, col="green")


lines(x = fpr, y = se, col='green')

lines(x = pfpr, y = pse, col='green')

#Analysis25
df = read.csv("F:/PHD/Projects/Stata/Madareg/Data/analysis25.csv")
df$se = with(df, tp/(tp + fn))
df$sp = with(df, tn/(tn + fp))
df$fpr = with(df, 1 - sp)
Sigma = matrix(c(.8091411, -.02225418, -.02225418, .20426229), ncol=2, byrow = FALSE)
Center = matrix(c(1.6773485, 2.314787, 2.3861435, 3.4766417, .82893074,
                  1.1267936, 1.3128566, 1.9194634, 2.1310296, 2.2901141,
                  2.6389507, 2.6296212), ncol=2, byrow=FALSE)
Vdf  = read.csv("F:/PHD/Projects/Stata/Madareg/Data/Vanalysis25.csv")
V = as.matrix(Vdf, ncol=12)

location = c("Africa", "Asia", "Central and South America", "Europe", "North America", "Oceania and Pacific")

for (i in 1:6){
  myplot(V=V[c(i,i+6), c(i, i+6)],
         Sigma = Sigma,
         Center = Center[i,],
         df = df[df$location==location[i],])
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                       CONFIDENCE INTERVALS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
df$se = with(df, TP/(TP + FN))
df$Dis = with(df, TP + FN)
df$seLower = NULL
df$seUpper = NULL
alpha = 0.05
for (r in 1:nrow(df)){
  df$seLower = ifelse(df$TP==0, 0, 1/(1 + (df$Dis - df$TP + 1)/(df$TP * qf(alpha/2, 2*df$TP, 2*(df$Dis - df$TP + 1)))))
  df$seUpper = ifelse(df$TP==df$Dis, 1, 1/(1 + (df$Dis - df$TP)/((df$TP + 1) * qf(1 - alpha/2, 2*(df$TP + 1), 2*(df$Dis - df$TP)))))
  
  
}

#=================================================
library(Metatron)


(telo<-fit.bivar(TP=TP,FN=FN,TN=TN,FP=FP,study=Study,data=df))
summary(telo)

(dementia<-fit.bivar(TP=TP,FN=FN,TN=TN,FP=FP,study=Study,data=df))
summary(dementia)


(ascus<-fit.bivar(TP=TP,FN=FN,TN=TN,FP=FP,mods=Test, study=StudyID,data=df))
summary(ascus)
