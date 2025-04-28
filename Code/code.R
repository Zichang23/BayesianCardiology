#read in data
df3 <- read.csv(file="Texas.csv") 
df3 <- df3[,-1]
head(df3)

############################
#     Data Exploration    
############################
#summary of dataset
summary(df3)

#pairs plot
library(GGally) 
ggpairs(df3, columns=5:7)

#boxplots of SMR overtime
df3$Year <- as.factor(df3$Year)
boxplot(SMR~Year, data=df3, ylab= "SMR", col="salmon")

#boxplots for obesity
boxplot(Obesity~Year, data=df3, ylab= "Obesity Rate (%)", col="lightblue")

#boxplots for diabetes
boxplot(Diabetes~Year, data=df3, ylab= "Diabetes Rate (%)", col="pink")

#boxplots for smoking
boxplot(Smoke~Year, data=df3, ylab= "Smoking Rate (%)", col="lightgoldenrod1")

#spatial map of avg SMR from 2011 to 2020
library(rgdal)
TX<- readOGR(dsn = "TX_210.shp")
by_TX<- group_by(df3, Code)
avgSMR <- summarize(by_TX, SMR = mean(SMR, na.rm=T))
avgSMR.TX<- merge(x=TX, y=avgSMR, by.x="GEOID", by.y="Code", all.x=FALSE)
library(leaflet)
avgSMR.TX.sp <- spTransform(avgSMR.TX, CRS("+proj=longlat +datum=WGS84 +no_defs")) 
myvar <- avgSMR.TX.sp@data$SMR
colors <- colorNumeric(palette = "YlOrBr", domain = myvar, reverse=FALSE)
leaflet(data=avgSMR.TX.sp) %>%
  addTiles() %>%
  addPolygons(fillColor =  ~colors(myvar), color="", fillOpacity = 0.8,
              weight = 1, smoothFactor = 0.55, opacity = 1.0) %>%
  addLegend(pal = colors, values = myvar, opacity = 1, title="SMR") %>% addScaleBar(position="bottomleft")

#spatial map with NA
library(rgdal)
TX<- readOGR(dsn = "TX_Counties.shp")
by_TX<- group_by(df3, Code)
avgSMR <- summarize(by_TX, SMR = mean(SMR, na.rm=T))
avgSMR.TX<- merge(x=TX, y=avgSMR, by.x="GEOID", by.y="Code", all.x=T)
library(leaflet)
avgSMR.TX.sp <- spTransform(avgSMR.TX, CRS("+proj=longlat +datum=WGS84 +no_defs")) 
myvar <- avgSMR.TX.sp@data$SMR
colors <- colorNumeric(palette = "YlOrBr", domain = myvar, reverse=FALSE)
leaflet(data=avgSMR.TX.sp) %>%
  addTiles() %>%
  addPolygons(fillColor =  ~colors(myvar), color="", fillOpacity = 0.8,
              weight = 1, smoothFactor = 0.55, opacity = 1.0) %>%
  addLegend(pal = colors, values = myvar, opacity = 1, title="SMR") %>% addScaleBar(position="bottomleft")

############################
#     Model Fitting    
############################
#model 1
#sample size
n <- nrow(df3)
#get values for independent variables
X <- df3[,c(6:8)]
#get values for coefficients and response variable y
myglm <- glm(formula = df3$Deaths ~ offset(log(df3$Expect)) + df3$Smoke + df3$Diabetes + df3$Obesity, family = "poisson")
beta <- c(myglm$coefficients[1], myglm$coefficients[2], myglm$coefficients[3], myglm$coefficients[4])
theta <- exp(beta[1] + beta[2] * X[,3] + beta[3]*X[,2] + beta[4]*X[,1])
set.seed(1) 
y <- rpois(n, lambda = theta * df3$Expect)

#likelihood
lik <- function(x){
  #set values for beta0, beta1, beta2, and beta3
  b0 <- x[1] 
  b1 <- x[2] 
  b2 <- x[3] 
  b3 <- x[4]
  theta <- exp(beta[1] + beta[2] * X[,3] + beta[3]*X[,2] + beta[4]*X[,1])
  # Loglikelihood function
  loglik <- sum(dpois(y, lambda = theta*df3$Expect, log=T)) 
  return(loglik)
}

#prior
log.prior <- function(x){
  #set values for beta0, beta1, beta2, and beta3
  b0 <- x[1] 
  b1 <- x[2] 
  b2 <- x[3] 
  b3 <- x[4]
  #calculate priors for beta0, beta1, beta2, and beta3
  b0.prior <- dnorm(b0, 0, 1, log=TRUE)
  b1.prior <- dnorm(b1, 0, 1, log=TRUE)
  b2.prior <- dnorm(b2, 0, 1, log=TRUE)
  b3.prior <- dnorm(b3, 0, 1, log=TRUE)
  #return the log priors
  return(b0.prior + b1.prior + b2.prior + b3.prior) 
}

#posterior
post <- function(x){
  return (lik(x) + log.prior(x)) 
}

#proposal
proposal <- function(x){
  return(rnorm(4, mean = x, sd = 0.1))
}

#Metropolis-Hastings sampling
myMCMC <- function(initial, iter){
  # create array to save accepted proposed values
  chain <- array(dim=c(iter + 1, 4)) 
  # set initial value
  chain[1, ] <- initial 
  for (i in 1:iter){
    # draw y from the proposal function
    y <- proposal(chain[i, ]) 
    # calculate the probability to accept the proposed values
    prob <- exp(post(y) - post(chain[i, ]))
    # determine whether to accept or reject Y 
    if (runif(1) < prob) {
      chain[i+1, ] <- y
    }else{ 
      chain[i+1, ] <- chain[i, ]
    }
  }
  return(chain)
}
#set initial values
initial <- beta
#set number of iterations
iter <- 22 * 10000

#generate MCMC chains
set.seed(1)
chain1 <- myMCMC(initial, iter)
chain2 <- myMCMC(initial, iter)
chain3 <- myMCMC(initial, iter)
#set thin=10
thin1 <- matrix(NA, ncol = 4, nrow = 10000)
for (i in 1:10000){
  if (i == 1){
    thin1[i, ] <- chain1[i*22,]
  } else {
    thin1[i, ] <- chain1[i*22,]
  }
}
thin2 <- matrix(NA, ncol = 4, nrow = 10000)
for (i in 1:10000){
  if (i == 1){
    thin2[i, ] <- chain2[i*22,]
  } else {
    thin2[i, ] <- chain2[i*22,]
  }
}
thin3 <- matrix(NA, ncol = 4, nrow = 10000)
for (i in 1:10000){
  if (i == 1){
    thin3[i, ] <- chain3[i*22,]
  } else {
    thin3[i, ] <- chain3[i*22,]
  }
}

#burn in the first 5000 of the iteration for each chain
burn <- 5000
c1 <- thin1[-c(1:burn), ]
c2 <- thin2[-c(1:burn), ]
c3 <- thin3[-c(1:burn), ]

#Check convergence
library(coda)
post.list=list(c1,c2,c3)
mcmc_list <- lapply(post.list, mcmc)
mcmc_list <- mcmc.list(mcmc_list)
gelman.diag(mcmc_list)
gelman.plot(mcmc_list)
plot(mcmc_list)

#calculate DIC for model 1
#method 1
beta.post <- rbind(c1,c2,c3)
aa <- colMeans(beta.post)
theta.pp <- exp(aa[1] + aa[2] * X[,3] + aa[3]*X[,2] + aa[4]*X[,1])
theta.mm <- mean(theta.pp)
logLike <- function(theta) sum(dpois(df3$Deaths, theta, log=T))
pDIC <- 2*(logLike(theta.mm) - mean(sapply(theta.pp, logLike) ))
-2*logLike(theta.mm) + 2*pDIC#42655.75
#method 2
library(CARBayes)
set.seed(1)
model1 <- S.glm(formula=Deaths ~offset(log(Expect)) + Smoke + Diabetes + Obesity, 
                family="poisson", data=df3, burnin=2000, n.sample=22000, thin=4)
model1$modelfit[1]#DIC = 47728.99

#determine spacial autocorrelation
myglm <- glm(formula=Deaths ~offset(log(Expect)) + Smoke + Diabetes + Obesity, family="poisson", data=df3) 
df3$res <- residuals(myglm) 
library("dplyr")
res2011 <- filter(df3, Year==2011) 
res2011.TX <- merge(x=TX, y=res2011, by.x="GEOID", by.y="Code", all.x=FALSE)

#calculate neighborhood matrix W
library(spdep)
W.poly <- poly2nb(res2011.TX, row.names = res2011.TX@data$GEOID) 
W <- nb2mat(W.poly, style = "B")
W.list <- nb2listw(W.poly, style = "B")

#calculate Moran’s I statistic
moran.mc(x = res2011.TX$res, listw = W.list, nsim = 10000, zero.policy=TRUE)

#determine temporal autocorrelation
pacf(df3$Deaths, main="Partial Autocorrelation Function")

#model 2
#order data
sp.order <- data.frame(Code=res2011.TX@data$GEOID, spatialorder=1:nrow(res2011.TX@data)) 
sp.new <- merge(x=df3, y=sp.order, by="Code") 
df4 <- arrange(sp.new, Year, spatialorder)
#fit model and generate three chains MCMC samples
library(CARBayesST)
chain1 <- ST.CARar(formula=Deaths ~offset(log(Expect)) + Smoke + Diabetes + Obesity, family="poisson",
                   data=df4, W=W, prior.var.beta=c(1,1,1,1), burnin=20000, n.sample=220000, thin=40, verbose=FALSE, AR=2)  
chain2 <- ST.CARar(formula=Deaths ~offset(log(Expect)) + Smoke + Diabetes + Obesity, family="poisson",
                   data=df4, W=W, prior.var.beta=c(1,1,1,1), burnin=20000, n.sample=220000, thin=40, verbose=FALSE, AR=2) 
chain3 <- ST.CARar(formula=Deaths ~offset(log(Expect)) + Smoke + Diabetes + Obesity, family="poisson",
                   data=df4, W=W, prior.var.beta=c(1,1,1,1), burnin=20000, n.sample=220000, thin=40, verbose=FALSE, AR=2)

#check convergence
library(coda)
beta.samples <- mcmc.list(chain1$samples$beta, chain2$samples$beta, chain3$samples$beta) 
plot(beta.samples)

#diagnostic
gelman.diag(beta.samples)
gelman.plot(beta.samples)

#summary of posterior results
print(chain1)

#check autocorrelation
autocorr.plot(beta.samples)

############################
#        Inference    
############################

#calculate relative risks
beta.samples.combined <- rbind(chain1$samples$beta[ ,2:4], chain2$samples$beta[ ,2:4], chain3$samples$beta[ ,2:4])
sd(df3$Smoke)
sd(df3$Diabetes)
sd(df3$Obesity)
#relative risk for obesity
round(quantile(exp(sd(df4$Obesity) * beta.samples.combined[ ,1]), c(0.5, 0.025, 0.975)),3)
#relative risk for diabetes
round(quantile(exp(sd(df4$Diabetes) * beta.samples.combined[ ,2]), c(0.5, 0.025, 0.975)),3)
#relative risk for smoking
round(quantile(exp(sd(df4$Smoke) * beta.samples.combined[ ,3]), c(0.5, 0.025, 0.975)),3)

#spatio-temporal trends in disease risk
#sample of fitted values
fit.post <- rbind(chain1$samples$fitted, chain2$samples$fitted, chain3$samples$fitted)
#sample size
n.samp <- nrow(fit.post)
#number of variables
n.var <- ncol(fit.post) 
#calculate posterior sample for relative risk
risk.post <- fit.post /matrix(rep(df4$Expect, n.samp), nrow=n.samp, ncol=n.var, byrow=TRUE)
#number of years
N <- length(table(df4$Year)) 
#posterior distribution of average spatial trends in relative risk
risk.change <- array(NA, c(n.samp, N))
for(i in 1:n.samp)
{
  risk.change[i, ] <- tapply(risk.post[i, ], df4$Year, mean) }

#calculate the 95% credible intervals for average spatially risk
trend <- as.data.frame(t(apply(risk.change, 2, quantile, c(0.5, 0.025, 0.975)))) 
trend <- trend %>% mutate(Year=names(table(df4$Year)))  
colnames(trend)[1:3] <- c("Median","LCI", "UCI")

#draw a plot to show the estimated temporal trend in disease risk
plot(trend$Year, trend$Median, xlab="Year", ylab="Median", type= "l", col="lightcoral",lwd=3)
lines(trend$Year, trend$LCI, lty=2, lwd=2, col="lightcoral")
lines(trend$Year, trend$UCI, lty=2, lwd=2, col="lightcoral")

#calculate the median risk and the posterior probability that the disease risks is greater than 1
risk2011 <- risk.post[ ,df4$Year==2011]  
risk2011.median <- apply(risk2011, 2, median) 
risk2011.high <- apply(risk2011 > 1, 2, mean) 
res2011.TX$risk2011.median <- c(risk2011.median)
res2011.TX$risk2011.high <- c(risk2011.high) 

#map for median risk
library(leaflet)
risk2011.TX.sp <- spTransform(res2011.TX, CRS("+proj=longlat +datum=WGS84 +no_defs"))
myvar2 <- risk2011.TX.sp@data$risk2011.median
colors2 <- colorNumeric(palette = "YlOrBr", domain = myvar2, reverse=FALSE)
leaflet(data=risk2011.TX.sp) %>%
  addTiles() %>%
  addPolygons(fillColor =  ~colors2(myvar2), color="", fillOpacity = 0.75,
              weight = 1, smoothFactor = 0.5, opacity = 1.0) %>%
  addLegend(pal = colors2, values = myvar2, opacity = 1, title="SMR") %>% addScaleBar(position="bottomleft")

#map for probability that the risk greater than one
myvar3 <- risk2011.TX.sp@data$risk2011.high
colors3 <- colorNumeric(palette = "YlOrBr", domain = myvar3, reverse=FALSE)
leaflet(data=risk2011.TX.sp) %>%
  addTiles() %>%
  addPolygons(fillColor =  ~colors3(myvar3), color="", fillOpacity = 0.9,
              weight = 1, smoothFactor = 0.5, opacity = 1.0) %>%
  addLegend(pal = colors3, values = myvar3, opacity = 1, title="SMR") %>% addScaleBar(position="bottomleft")

#health inequalities
risk.median <- apply(risk.post, 2, median)
unequal <- tapply(risk.median, df4$Year, IQR)

#draw the health inequalities trend over years
plot(trend$Year, inequality, xlab="Year", ylab="Inequality", type= "l", col="lightcoral",lwd=3)
