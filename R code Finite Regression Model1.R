library(dplyr)
library(xtable)
library(fmrs)
library(nnet)
library(ggplot2)
library(reshape2)
####Simulation####

set.seed(1980)
nComp = 2
nCov = 10
nObs = 1000
dispersion = c(1, 1)
mixProp = c(0.3, 0.7)
rho = 0.5
coeff1 = c( 1,  0, 0, 3, 0, 1, 0.6, 0,  3, 0,  1)
coeff2 = c(-1, 2,  0,  0, 3, -1, 0, 0, 4, .7, -2)
umax = 40
dat <- fmrs.gendata(nObs = nObs, nComp = nComp, nCov = nCov,
                    coeff = c(coeff1, coeff2), dispersion = dispersion,
                    mixProp = mixProp, rho = rho, umax = umax,
                    disFamily = "norm")

res.mle <- fmrs.mle(y = dat$y, x = dat$x, delta = dat$delta,
                    nComp = nComp, disFamily = "norm",
                    initCoeff = rnorm(nComp*nCov+nComp),
                    initDispersion = rep(1, nComp),
                    initmixProp = rep(1/nComp, nComp))
res.mle
x1<-coefficients(res.mle)
xtable(x1,digits=3,caption="Maximum Likelihood Estimate of Finite Mixture of Regression Model")
###Mixlasso###
res.lam1 <- fmrs.tunsel(y = dat$y, x = dat$x, delta = dat$delta,
                        nComp = nComp, disFamily = "norm",
                        initCoeff = c(coefficients(res.mle)),
                        initDispersion = dispersion(res.mle),
                        initmixProp = mixProp(res.mle),
                        penFamily = "lasso")
res.var1 <- fmrs.varsel(y = dat$y, x = dat$x, delta = dat$delta,
                        nComp = ncomp(res.mle), disFamily = "norm",
                        initCoeff=c(coefficients(res.mle)),
                        initDispersion = dispersion(res.mle),
                        initmixProp = mixProp(res.mle),
                        penFamily = "lasso",
                        lambPen = slot(res.lam1, "lambPen"))



x2<-round(coefficients(res.var1)[-1,],5)


###mixScad####

res.lam2 <- fmrs.tunsel(y = dat$y, x = dat$x, delta = dat$delta,
                        nComp = nComp, disFamily = "norm",
                        initCoeff = c(coefficients(res.mle)),
                        initDispersion = dispersion(res.mle),
                        initmixProp = mixProp(res.mle),
                        penFamily = "scad")
res.var2 <- fmrs.varsel(y = dat$y, x = dat$x, delta = dat$delta,
                        nComp = ncomp(res.mle), disFamily = "norm",
                        initCoeff=c(coefficients(res.mle)),
                        initDispersion = dispersion(res.mle),
                        initmixProp = mixProp(res.mle),
                        penFamily = "scad",
                        lambPen = slot(res.lam2, "lambPen"))



x3<-round(coefficients(res.var2)[-1,],5)


###Mixhard####
res.lam3 <- fmrs.tunsel(y = dat$y, x = dat$x, delta = dat$delta,
                        nComp = nComp, disFamily = "norm",
                        initCoeff = c(coefficients(res.mle)),
                        initDispersion = dispersion(res.mle),
                        initmixProp = mixProp(res.mle),
                        penFamily = "hard")
res.var3 <- fmrs.varsel(y = dat$y, x = dat$x, delta = dat$delta,
                        nComp = ncomp(res.mle), disFamily = "norm",
                        initCoeff=c(coefficients(res.mle)),
                        initDispersion = dispersion(res.mle),
                        initmixProp = mixProp(res.mle),
                        penFamily = "hard",
                        lambPen = slot(res.lam3, "lambPen"))



x4<-round(coefficients(res.var3)[-1,],5)
xtable(x2,digits=3,caption="Penalized Maximum Likelihood Estimate of Finite Mixture of Regression Model(MIXLASSO)")

xtable(x3,digits=3,caption="Penalized Maximum Likelihood Estimate of Finite Mixture of Regression Model(MIXSCAD)")

xtable(x4,digits=3,caption="Penalized Maximum Likelihood Estimate of Finite Mixture of Regression Model(MIXHARD)")
dispersion(res.var1)
mixProp(res.var1)
dispersion(res.var2)
mixProp(res.var2)
dispersion(res.var3)
mixProp(res.var3)


#### Breast cancer dataset####

ml<-read.csv("breastc.csv")
head(ml)
z<-as.factor(ml$diagnosis)
str(z)
data_y=as.numeric(z)

data_x<- model.matrix( ~ ml$radius_mean+ml$texture_mean+ml$perimeter_mean+ml$area_mean+ml$smoothness_mean+ml$compactness_mean+
                         ml$concavity_mean+ml$concave.points_mean+
                         ml$symmetry_mean+ml$fractal_dimension_mean+ml$radius_se+ml$texture_se+
                         ml$perimeter_se+ml$area_se+ml$smoothness_se+ml$compactness_se+ml$concavity_se+
                         ml$concave.points_se+ml$symmetry_se+ml$fractal_dimension_se+ml$radius_worst+ml$texture_worst+
                         ml$perimeter_worst+ml$area_worst+ml$smoothness_worst+ml$compactness_worst+ml$concavity_worst+
                         ml$concave.points_worst+ml$symmetry_worst+ml$fractal_dimension_worst)[,-1]

res.mle <- fmrs.mle(y = data_y, x = data_x, delta = ml$delta,
                    nComp = 2, disFamily = "norm",
                    initCoeff = rnorm(nComp*2+nComp),
                    initDispersion = rep(1, nComp),
                    initmixProp = rep(1/nComp, nComp))

res.lam1 <- fmrs.tunsel(y = data_y, x = data_x, delta = ml$delta,
                        nComp = 2, disFamily = "norm",
                        initCoeff = c(coefficients(res.mle)),
                        initDispersion = dispersion(res.mle),
                        initmixProp = mixProp(res.mle),
                        penFamily = "lasso")
res.var1 <- fmrs.varsel(y = data_y, x = data_x, delta = data$delta,
                        nComp = ncomp(res.mle), disFamily = "norm",
                        initCoeff=c(coefficients(res.mle)),
                        initDispersion = dispersion(res.mle),
                        initmixProp = mixProp(res.mle),
                        penFamily = "lasso",
                        lambPen = slot(res.lam1, "lambPen"))


x1<-coefficients(res.mle)
x2<-round(coefficients(res.var1)[-1,],5)
xtable(x1,digits=3,caption="Maximum Likelihood Estimate of Finite Mixture of Regression Model (Breast Cancer Data)")
xtable(x2,digits=3,caption="Maximum Likelihood Estimate of Finite Mixture of Regression Model(MIXLASSO)")
###MIXSCAD###

res.lam2 <- fmrs.tunsel(y = data_y, x = data_x, delta = ml$delta,
                        nComp = 2, disFamily = "norm",
                        initCoeff = c(coefficients(res.mle)),
                        initDispersion = dispersion(res.mle),
                        initmixProp = mixProp(res.mle),
                        penFamily = "scad")
res.var2 <- fmrs.varsel(y = data_y, x = data_x, delta = data$delta,
                        nComp = ncomp(res.mle), disFamily = "norm",
                        initCoeff=c(coefficients(res.mle)),
                        initDispersion = dispersion(res.mle),
                        initmixProp = mixProp(res.mle),
                        penFamily = "scad",
                        lambPen = slot(res.lam2, "lambPen"))
x3<-round(coefficients(res.var2)[-1,],5)
xtable(x3,digits=3,caption="Maximum Likelihood Estimate of Finite Mixture of Regression Model(MIXSCAD)")
###MIXHARD###
res.lam3 <- fmrs.tunsel(y = data_y, x = data_x, delta = ml$delta,
                        nComp = 2, disFamily = "norm",
                        initCoeff = c(coefficients(res.mle)),
                        initDispersion = dispersion(res.mle),
                        initmixProp = mixProp(res.mle),
                        penFamily = "scad")
res.var3 <- fmrs.varsel(y = data_y, x = data_x, delta = data$delta,
                        nComp = ncomp(res.mle), disFamily = "norm",
                        initCoeff=c(coefficients(res.mle)),
                        initDispersion = dispersion(res.mle),
                        initmixProp = mixProp(res.mle),
                        penFamily = "scad",
                        lambPen = slot(res.lam3, "lambPen"))
x4<-round(coefficients(res.var3)[-1,],5)
xtable(x3,digits=3,caption="Maximum Likelihood Estimate of Finite Mixture of Regression Model(MIXHARD)")



#### heart data set#####
df<-read.csv("heart.csv")
head(ml)
df$sex <- as.factor(df$sex)
df$cp <- as.factor(df$cp)
df$fbs <- as.factor(df$fbs)
df$restecg <- as.factor(df$restecg)
df$exang <- as.factor(df$exang)
df$slope <- as.factor(df$slope)
df$ca <- as.factor(df$ca)
df$thal <- as.factor(df$thal)
df$target <- as.factor(df$target)
y1<-as.numeric(y)
x <- model.matrix(target~., df)[,-1]


res.mle <- fmrs.mle(y = y1, x = x, delta = df$delta,
                    nComp = 2, disFamily = "norm",
                    initCoeff = rnorm(nComp*2+nComp),
                    initDispersion = rep(1, nComp),
                    initmixProp = rep(1/nComp, nComp))
coefficients(res.mle)
res.lam1 <- fmrs.tunsel(y = y1, x = x, delta = df$delta,
                        nComp = 2, disFamily = "norm",
                        initCoeff = c(coefficients(res.mle)),
                        initDispersion = dispersion(res.mle),
                        initmixProp = mixProp(res.mle),
                        penFamily = "lasso")
res.var1 <- fmrs.varsel(y = y1, x = x, delta = df$delta,
                        nComp = ncomp(res.mle), disFamily = "norm",
                        initCoeff=c(coefficients(res.mle)),
                        initDispersion = dispersion(res.mle),
                        initmixProp = mixProp(res.mle),
                        penFamily = "lasso",
                        lambPen = slot(res.lam1, "lambPen"))


x1<-coefficients(res.mle)
x2<-round(coefficients(res.var1)[-1,],5)
xtable(x1,digits=3,caption="Maximum Likelihood Estimate of Finite Mixture of Regression Model (Heart data)")

###MIXSCAD###

res.lam2 <- fmrs.tunsel(y = y1, x = x, delta = df$delta,
                        nComp = 2, disFamily = "norm",
                        initCoeff = c(coefficients(res.mle)),
                        initDispersion = dispersion(res.mle),
                        initmixProp = mixProp(res.mle),
                        penFamily = "scad")
res.var2 <- fmrs.varsel(y =y1, x = x, delta = df$delta,
                        nComp = ncomp(res.mle), disFamily = "norm",
                        initCoeff=c(coefficients(res.mle)),
                        initDispersion = dispersion(res.mle),
                        initmixProp = mixProp(res.mle),
                        penFamily = "scad",
                        lambPen = slot(res.lam2, "lambPen"))
x3<-round(coefficients(res.var2)[-1,],5)

###MIXHARD###
res.lam3 <- fmrs.tunsel(y = y1, x = x, delta = df$delta,
                        nComp = 2, disFamily = "norm",
                        initCoeff = c(coefficients(res.mle)),
                        initDispersion = dispersion(res.mle),
                        initmixProp = mixProp(res.mle),
                        penFamily = "scad")
res.var3 <- fmrs.varsel(y = y1, x = x, delta = df$delta,
                        nComp = ncomp(res.mle), disFamily = "norm",
                        initCoeff=c(coefficients(res.mle)),
                        initDispersion = dispersion(res.mle),
                        initmixProp = mixProp(res.mle),
                        penFamily = "scad",
                        lambPen = slot(res.lam3, "lambPen"))
x4<-round(coefficients(res.var3)[-1,],5)
xtable(x2,digits=3,caption="Maximum Likelihood Estimate of Finite Mixture of Regression Model(MIXLASSO)")
xtable(x3,digits=3,caption="Maximum Likelihood Estimate of Finite Mixture of Regression Model(MIXSCAD)")
xtable(x3,digits=3,caption="Maximum Likelihood Estimate of Finite Mixture of Regression Model(MIXHARD)")

