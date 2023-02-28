
#Data
ying1 <- read.csv('https://raw.githubusercontent.com/VNyaga/ESMARConf2023/main/Data/Ying1.csv', sep=',')

#metafor
#============================================================================================
library(metafor)
ying1$dis = with(ying1, tp + fn)
ying1$ndis = with(ying1, tn + fp)

#fit logistic-normal model
se = rma.glmm(xi=tp, ni=dis, data=ying1, measure ="PLO")
sp = rma.glmm(xi=tn, ni=ndis, data=ying1, measure ="PLO")

#Results
sp
forest.rma(sp)

#Transform to 0-1 scale
c(exp(sp$b)/(1 + exp(sp$b)), exp(sp$ci.lb)/(1 + exp(sp$ci.lb)), exp(sp$ci.ub)/(1 + exp(sp$ci.ub)))

#fit logistic model
se = rma.glmm(xi=tp, ni=dis, data=ying1, measure ="PLO", method="EE")

#Results
se
forest.rma(se)

#Transform to 0-1 scale
c(exp(se$b)/(1 + exp(se$b)), exp(se$ci.lb)/(1 + exp(se$ci.lb)), exp(se$ci.ub)/(1 + exp(se$ci.ub)))


#MADA
#============================================================================================
library(mada)

#names should be in capitals
library(stringr)
names(ying1)[2:5] <- str_to_upper(names(ying1)[2:5])

#fit normal-normal model
fit = reitsma(ying1, correction = 0)

#fit normal-normal model with default 0.5 continuity correction
fit = reitsma(ying1)

#forest plot
par(mfcol=c(1,2))
forest(madad(ying1), type = "sens", snames=ying1$author, xlab="Sens")
forest(madad(ying1), type = "spec", snames=ying1$author, xlab="Spec")

#Transform to 0-1 scale
ss = summary(fit)

c(ss$coefficients[3,1], ss$coefficients[3,5], ss$coefficients[3,6]) #Sens
c(1-ss$coefficients[4,1], 1-ss$coefficients[4,6], 1-ss$coefficients[4,5]) #Spec


#CopulaDTA
#============================================================================================
library(CopulaDTA)

#Studyid should be a factor variable
ying1$author = as.factor(ying1$author)

#Define the models i.e. choose the copula
gauss <-  cdtamodel(copula = 'gauss') #models full-range correlation
gauss #see the model specification

fgm <-  cdtamodel(copula = 'fgm') #models weak correlation
fgm #see the model specification

#Fit GAUSSIAN copula with normal marginals
T0 = Sys.time()
fit1 <- fit(gauss,
            SID='author',
            data=ying1,
            iter=2000,
            warmup=1000,
            thin=1,
            seed=3)

T1 = Sys.time() - T0; T1
traceplot(fit1)

print(fit1)

plot(fit1)

#Fit FGM copula with beta marginals
T0 = Sys.time()
fit2 <- fit(fgm,
            SID='author',
            data=ying1,
            iter=2000,
            warmup=1000,
            thin=1,
            seed=3)

T1 = Sys.time() - T0; T1

traceplot(fit2)

print(fit2)

plot(fit2)

