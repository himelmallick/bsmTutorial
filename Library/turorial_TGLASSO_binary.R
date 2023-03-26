###################
# Clear workspace #
###################

rm(list = ls())

#########################
# Set working directory #
#########################

workingDirectory = getwd()
setwd(workingDirectory)

##################
# Load libraries #
##################

library(tidyverse) # for data manipulation
library(bayesreg) # For fitting the model
library(plotrix) # For CI plot
library(coda) # for MCMC data
library(pROC) # For AUC calculation
library(Rfssa) # To load github data

#############
# Load data #
#############

load_github_data("https://github.com/himelmallick/bspTutorial/blob/master/Data/TGLASSO_filtered.RData")
trainX = pcl$trainX
trainY = pcl$trainY

#########################
# Remove missing values #
#########################

missing_index = which(is.na(trainY$tamoxifen))
trainY = trainY[-missing_index,]
trainX = trainX[-missing_index, ]

################
# Sanity Check #
################

all(rownames(trainX) == rownames(trainY))

###################
# Standardization #
###################

y.train = trainY$tamoxifen - mean(trainY$tamoxifen)
x.train = scale(as.matrix(trainX))
df.train = cbind.data.frame(y.train, x.train)

#####################
# Creating binary Y #
#####################

y.train.binary = ifelse(trainY$tamoxifen >= median(trainY$tamoxifen), 1, 0)
y.train.binary = as.factor(y.train.binary)
df.train=cbind.data.frame(y.train.binary, x.train)

#############
# Horseshoe #
#############

set.seed(1234)
fit.hs = bayesreg(y.train.binary ~., df.train, model="logistic", prior = "hs", burnin = 10000)

##################
# Horseshoe plus #
##################

fit.hsplus = bayesreg(y.train.binary ~., df.train, model="logistic", prior="horseshoe+", burnin = 10000)

##################
# Bayesian ridge #
##################

fit.br = bayesreg(y.train.binary ~., df.train, model="logistic", prior = "ridge", burnin = 10000)

##################
# Bayesian LASSO #
##################

fit.bl = bayesreg(y.train.binary ~., df.train, model="logistic", prior="lasso", burnin = 10000)

##########
# Result #
##########

################
# AUC/Accuracy #
################

#############
# Horseshoe #
#############

y.pred = predict(fit.hs, df.train, type='class')
table.val = table(y.train.binary, y.pred)
roc.obj = roc(as.numeric(y.train.binary), as.numeric(y.pred))
auc.val= auc(roc.obj)
accuracy.val = sum(diag(table.val))/sum(table.val)
measure.valhs  = c(auc.val, accuracy.val)

#################
# Horseshoeplus #
#################

y.pred = predict(fit.hsplus, df.train, type='class')
table.val = table(y.train.binary, y.pred)
roc.obj = roc(as.numeric(y.train.binary), as.numeric(y.pred))
auc.val= auc(roc.obj)
accuracy.val = sum(diag(table.val))/sum(table.val)
measure.valhsplus  = c(auc.val, accuracy.val)

##################
# Bayesian ridge #
##################

y.pred = predict(fit.br, df.train, type='class')
table.val = table(y.train.binary, y.pred)
roc.obj = roc(as.numeric(y.train.binary), as.numeric(y.pred))
auc.val= auc(roc.obj)
accuracy.val = sum(diag(table.val))/sum(table.val)
measure.valbr  = c(auc.val, accuracy.val)

##################
# Bayesian LASSO #
##################

y.pred = predict(fit.bl, df.train, type='class')
table.val = table(y.train.binary, y.pred)
roc.obj = roc(as.numeric(y.train.binary), as.numeric(y.pred))
auc.val= auc(roc.obj)
accuracy.val = sum(diag(table.val))/sum(table.val)
measure.valbl = c(auc.val, accuracy.val)

measure.val = rbind(measure.valhs[1], measure.valhsplus[1], measure.valbr[1],
                    measure.valbl[1])
colnames(measure.val) = c("AUC")
row.names(measure.val) = c("Horseshoe", "Horseshoe+", "Bayesian Ridge", "Bayesian LASSO")
measure.val


###################
# Save the output #
###################

save(fit.hs, file = paste(workingDirectory, './RdataOutput/fit.HSBinary', '.RData', sep = ''))
save(fit.hsplus, file = paste(workingDirectory, './RdataOutput/fit.HSplusBinary', '.RData', sep = ''))
save(fit.br, file = paste(workingDirectory, './RdataOutput/fit.BRBinary', '.RData', sep = ''))
save(fit.bl, file = paste(workingDirectory, './RdataOutput/fit.BLBinary', '.RData', sep = ''))


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

pdf("./Plots/CIplot_binary.pdf")
plotCI(c(m.cre1, m.cre2, m.cre3, m.cre4),
       x = c((1:Q)-2*xOff, (1:Q), (1:Q)+2*xOff, (1:Q)+4*xOff),
       ui = c(U1, U2, U3, U4),
       li = c(L1, L2, L3, L4),
       xlab = "betas",
       ylab = "Estimates",
       axes = FALSE,
       lwd = 1,
       ylim = c(-1, 1),
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
legend(8, 1, c("HS", "HS+", "BR", "BL"), pch = c(15, 16, 17, 18), pt.cex = 1,
       text.width = 1, col = c("violet", "brown", "purple", "deeppink"))
dev.off()

