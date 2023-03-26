###################
# Clear workspace #
###################

rm(list = ls())

##################
# Load libraries #
##################

library(tidyverse) # For data manipulation
library(bayesreg) # For fitting the model
library(plotrix) # For CI plot
library(coda) # for MCMC data
library(psych) # For histogram plot
library(survC1) # For calculation C-index for survival data
library(Rfssa) # To load github data

#########################
# Set working directory #
#########################

workingDirectory = getwd()
setwd(workingDirectory)

#############
# Load data #
#############

load_github_data("https://github.com/himelmallick/bspTutorial/blob/master/Data/Survival.RData")
trainX = pcl$trainX
trainY = c(pcl$trainY)

################
# Sanity Check #
################

all(rownames(trainX) == rownames(trainY))

###################
# Standardization #
###################

y.train = trainY$trainY - mean(trainY$trainY)
x.train = scale(as.matrix(trainX))
df.train = cbind.data.frame(y.train, x.train)

#############
# Horseshoe #
#############

set.seed(1234)
fit.hs = bayesreg(y.train ~., df.train, prior = "hs", burnin = 10000)
yhat = as.matrix(pcl$original_X_standardized)%*%as.matrix(fit.hs$beta)
mydata = data.frame(as.matrix(pcl$original_Y), yhat)
out = Est.Cval(mydata, tau = 2000, nofit = TRUE)
cindex.valhs = c(out$Dhat)

##############
# Horseshoe+ #
##############

fit.hsplus = bayesreg(y.train ~., df.train, prior = "hs+", burnin = 10000)
yhat = as.matrix(pcl$original_X_standardized)%*%as.matrix(fit.hsplus$beta)
mydata = data.frame(as.matrix(pcl$original_Y), yhat)
out = Est.Cval(mydata, tau = 2000, nofit = TRUE)
cindex.valhsplus = c(out$Dhat)

##################
# Bayesian ridge #
##################

fit.br = bayesreg(y.train ~., df.train, prior = "ridge", burnin = 10000)
yhat = as.matrix(pcl$original_X_standardized)%*%as.matrix(fit.br$beta)
mydata = data.frame(as.matrix(pcl$original_Y), yhat)
out = Est.Cval(mydata, tau = 2000, nofit = TRUE)
cindex.valbr =   c(out$Dhat)

##################
# Bayesian LASSO #
##################

fit.bl = bayesreg(y.train ~., df.train, prior = "lasso", burnin = 10000) # Lasso
yhat = as.matrix(pcl$original_X_standardized)%*%as.matrix(fit.bl$beta)
mydata = data.frame(as.matrix(pcl$original_Y), yhat)
out = Est.Cval(mydata, tau = 2000, nofit = TRUE)
cindex.valbl = c(out$Dhat)

##########
# Result #
##########

cindex.val = rbind(cindex.valhs, cindex.valhsplus, cindex.valbr, cindex.valbl)
colnames(cindex.val) = "C-Index"
row.names(cindex.val) = c("Horseshoe", "Horseshoe+", "Bayesian ridge", "Bayesian LASSO")
cindex.val

###################
# Save the output #
###################

save(fit.hs, file = paste(workingDirectory, './RdataOutput/fit.HSSurvival', '.RData', sep = ''))
save(fit.hsplus, file = paste(workingDirectory, './RdataOutput/fit.HSplusSurvival', '.RData', sep = ''))
save(fit.br, file = paste(workingDirectory, './RdataOutput/fit.BRSurvival', '.RData', sep = ''))
save(fit.bl, file = paste(workingDirectory, './RdataOutput/fit.BLSurvival', '.RData', sep = ''))


#################
# Visualization #
#################

###############
# Top 10 beta #
###############

# Find the absolute value of the posterior median
top_n = 10

# HS
topbeta = apply(abs(t(fit.hs$beta)), 2, median)
topbetas = topbeta[order(topbeta, decreasing = TRUE)]
index_top = rep(NA, top_n)
for(i in 1:top_n){
  index_top[i] = which(topbeta == topbetas[i])
}

# Extract top 10 beta
topBeta_hs = t(fit.hs$beta[index_top,])

# HS+
topbeta = apply(abs(t(fit.hsplus$beta)), 2, median)
topbetas = topbeta[order(topbeta, decreasing = TRUE)]
index_top = rep(NA, top_n)
for(i in 1:top_n){
  index_top[i] = which(topbeta == topbetas[i])
}

# Extract top 10 beta
topBeta_hsplus = t(fit.hsplus$beta[index_top,])


# BR
topbeta = apply(abs(t(fit.br$beta)), 2, median)
topbetas = topbeta[order(topbeta, decreasing = TRUE)]
index_top = rep(NA, top_n)
for(i in 1:top_n){
  index_top[i] = which(topbeta == topbetas[i])
}

# Extract top 10 beta
topBeta_br = t(fit.br$beta[index_top,])

# BL
topbeta = apply(abs(t(fit.bl$beta)), 2, median)
topbetas = topbeta[order(topbeta, decreasing = TRUE)]
index_top = rep(NA, top_n)
for(i in 1:top_n){
  index_top[i] = which(topbeta == topbetas[i])
}

# Extract top 10 beta
topBeta_bl = t(fit.bl$beta[index_top,])


###########
# CI plot #
###########

# Each posterior has p columns and M samples

cred1 = apply(topBeta_hs, 2, quantile, prob = c(0.025, 0.5, 0.975)) # HS
cred2 = apply(topBeta_hsplus, 2, quantile, prob = c(0.025, 0.5, 0.975)) # HS+
cred3 = apply(topBeta_br, 2, quantile, prob = c(0.025, 0.5, 0.975)) # BR
cred4 = apply(topBeta_bl, 2, quantile, prob = c(0.025, 0.5, 0.975)) # BL

L1 = cred1[1,]; L2 = cred2[1,]; L3 = cred3[1,]; L4 = cred4[1,]
m.cre1 = cred1[2,]; m.cre2 = cred2[2,]; m.cre3 = cred3[2,]; m.cre4 =cred4[2,]
U1 = cred1[3,]; U2 = cred2[3,]; U3 = cred3[3,]; U4 = cred4[3,]
xOff = 0.1
Q = 10

pdf("./Plots/CIplot_survival.pdf")
plotCI(c(m.cre1, m.cre2, m.cre3, m.cre4),
       x = c((1:Q)-2*xOff, (1:Q), (1:Q)+2*xOff, (1:Q)+4*xOff),
       ui = c(U1, U2, U3, U4),
       li = c(L1, L2, L3, L4),
       xlab = "betas",
       ylab = "Estimates",
       axes = FALSE,
       lwd = 1,
       ylim = c(-0.4, 0.45),
       font = 1,
       cex.lab = 1,
       scol = "orange3",
       col = rep(c("violet", "brown", "purple", "deeppink"), each = Q),
       cex = 1,
       cex.axis = 1,
       pch = rep(15:18, each = Q))
axis(1, 1:Q, font = 1, cex.axis = 1, adj = 0, las = 1)
axis(2)
box()
abline(h = 0, lty = 2, lwd = 2, col = "gray60")
legend(8, 0.45, c("HS", "HS+", "BR", "BL"), pch = c(15, 16, 17, 18), pt.cex = 1,
       text.width = 1, col = c("violet", "brown", "purple", "deeppink"))
dev.off()
