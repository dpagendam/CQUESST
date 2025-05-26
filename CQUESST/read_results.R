
library(cmdstanr)
library(rstan)
library(abind)
library(invgamma)

data1 <- read_cmdstan_csv(
  files = "output.csv",
  variables = c("D", "R", "F", "S", "H", "IOMCMInit", "sigma2_RPM", "sigma2_ROC", "sigma2_TOC", "DcmpRateDPM", "DcmpRateRPM", "DcmpRateBIOFv263", "DcmpRateBIOSv263", "DcmpRateHUM", "sigma2_D", "sigma2_R", "sigma2_F", "sigma2_S", "sigma2_H", "RatioDPMToRPM", "FracClay", "FracSolidXXXToBIOFv263", 
  "treatment_alpha_Fallow", "treatment_alpha_Pasture", "treatment_alpha_NN0", "treatment_alpha_NN1", "treatment_alpha_MM0", "treatment_alpha_MM1", "treatment_alpha_II0", "treatment_alpha_II1", "treatment_alpha_IM0", "treatment_alpha_IM1", "treatment_alpha_IN0", "treatment_alpha_IN1", "treatment_alpha_MN0", "treatment_alpha_MN1",
  "FracManureDPM", "FracManureRPM", "FracManureBIOF", "FracManureBIOS", "FracManureHUM", "FracSolidHUMToBIOSv263", "EvapFactor"),
  sampler_diagnostics = NULL,
  format = getOption("cmdstanr_draws_format", NULL)
)


draws1 <- data1$post_warmup_draws[seq(1, 10000, 10), , ]


sigma2_RPM_1 <- as.numeric(draws1[, , "sigma2_RPM"])

sigma2_ROC_1 <- as.numeric(draws1[, , "sigma2_ROC"])

sigma2_TOC_1 <- as.numeric(draws1[, , "sigma2_TOC"])

DcmpRateDPM_1 <- as.numeric(draws1[, , "DcmpRateDPM"])

DcmpRateRPM_1 <- as.numeric(draws1[, , "DcmpRateRPM"])

DcmpRateBIOFv263_1 <- as.numeric(draws1[, , "DcmpRateBIOFv263"])

DcmpRateBIOSv263_1 <- as.numeric(draws1[, , "DcmpRateBIOSv263"])

DcmpRateHUM_1 <- as.numeric(draws1[, , "DcmpRateHUM"])

sigma2_D_1 <- as.numeric(draws1[, , "sigma2_D"])

sigma2_R_1 <- as.numeric(draws1[, , "sigma2_R"])

sigma2_F_1 <- as.numeric(draws1[, , "sigma2_F"])

sigma2_S_1 <- as.numeric(draws1[, , "sigma2_S"])

sigma2_H_1 <- as.numeric(draws1[, , "sigma2_H"])

RatioDPMToRPM_1 <- as.numeric(draws1[, , "RatioDPMToRPM"])

FracClay_1 <- as.numeric(draws1[, , "FracClay"])

FracManureDPM_1 <- as.numeric(draws1[, , "FracManureDPM"])

FracManureRPM_1 <- as.numeric(draws1[, , "FracManureRPM"])

FracManureBIOF_1 <- as.numeric(draws1[, , "FracManureBIOF"])

FracManureBIOS_1 <- as.numeric(draws1[, , "FracManureBIOS"])

FracManureHUM_1 <- as.numeric(draws1[, , "FracManureHUM"])

FracSolidXXXToBIOFv263_1 <- as.numeric(draws1[, , "FracSolidXXXToBIOFv263"])

FracSolidHUMToBIOSv263_1 <- as.numeric(draws1[, , "FracSolidHUMToBIOSv263"])

treatment_alpha_Pasture_1 <- as.numeric(draws1[, , "treatment_alpha_Pasture"])

treatment_alpha_NN0_1 <- as.numeric(draws1[, , "treatment_alpha_NN0"])

treatment_alpha_NN1_1 <- as.numeric(draws1[, , "treatment_alpha_NN1"])

treatment_alpha_MM0_1 <- as.numeric(draws1[, , "treatment_alpha_MM0"])

treatment_alpha_MM1_1 <- as.numeric(draws1[, , "treatment_alpha_MM1"])

treatment_alpha_II0_1 <- as.numeric(draws1[, , "treatment_alpha_II0"]) 

treatment_alpha_II1_1 <- as.numeric(draws1[, , "treatment_alpha_II1"]) 

treatment_alpha_IM0_1 <- as.numeric(draws1[, , "treatment_alpha_IM0"]) 

treatment_alpha_IM1_1 <- as.numeric(draws1[, , "treatment_alpha_IM1"]) 

treatment_alpha_IN0_1 <- as.numeric(draws1[, , "treatment_alpha_IN0"])

treatment_alpha_IN1_1 <- as.numeric(draws1[, , "treatment_alpha_IN1"])
 
treatment_alpha_MN0_1 <- as.numeric(draws1[, , "treatment_alpha_MN0"])

treatment_alpha_MN1_1 <- as.numeric(draws1[, , "treatment_alpha_MN1"])

treatment_alpha_Fallow_1 <- as.numeric(draws1[, , "treatment_alpha_Fallow"])








# CroppingParameters (D): Range 5 - 20
# CroppingParameters (H): 0 - 0.5
# Cropping Parameters (R): 0 - 5

#TillageParameters: 0-1

constrainDensity <- function(xVals, yVals, lowerBound, upperBound)
{
	lowerInds <- which(xVals < lowerBound)
	upperInds <- which(xVals > upperBound)
	yVals[c(lowerInds, upperInds)] <- 0
	return(yVals)
}



#Priors and Posteriors on Decay Rates


pdf("DecayRate_D.pdf", height = 4, width = 4)
goodInds <- which(!is.na(DcmpRateDPM_1))
K_D <- density(DcmpRateDPM_1[goodInds])
x <- c(K_D$x)
y <- c(K_D$y)
plot(0, 0, col = "white", xlim = c(5, 20), ylim = c(0, 1.2*max(y)), xlab = "Decay Rate", ylab = "Density", main = "D Pool")
x = c(K_D$x, rev(K_D$x), K_D$x[1])
y = c(rep(0, length(K_D$y)), rev(K_D$y), 0)
y <- constrainDensity(x, y, 5, 20)
polygon(x, y, col = adjustcolor("purple", 0.25), border = "purple")

x_seq <- seq(5, 20, by = 0.01)
y_seq <- dnorm(x_seq, 10.0, 0.5)/(1 - pnorm(5, 10, 0.5))

prior_x <- c(x_seq, rev(x_seq), x_seq[1])
prior_y <- c(y_seq, rep(0, length(y_seq)), y_seq[1])
polygon(prior_x, prior_y, col = adjustcolor("grey", 0.25), border = "grey")
dev.off()


pdf("DecayRate_R.pdf", height = 4, width = 4)
goodInds <- which(!is.na(DcmpRateRPM_1))
K_R <- density(DcmpRateRPM_1[goodInds])
x <- c(K_R$x)
y <- c(K_R$y)
plot(0, 0, col = "white", xlim = c(0.05, 0.15), ylim = c(0, 1.2*max(y)), xlab = "Decay Rate", ylab = "Density", main = "R Pool")
x = c(K_R$x, rev(K_R$x), K_R$x[1])
y = c(rep(0, length(K_R$y)), rev(K_R$y), 0)
y <- constrainDensity(x, y, 0.05, 5)
polygon(x, y, col = adjustcolor("purple", 0.25), border = "purple")

x_seq <- seq(0.05, 5, by = 0.0001)
y_seq <- dnorm(x_seq, 0.07, 0.0035)/(1 - pnorm(0.05, 0.07, 0.0035))

prior_x <- c(x_seq, rev(x_seq), x_seq[1])
prior_y <- c(y_seq, rep(0, length(y_seq)), y_seq[1])
polygon(prior_x, prior_y, col = adjustcolor("grey", 0.25), border = "grey")
dev.off()




pdf("DecayRate_F.pdf", height = 4, width = 4)
goodInds <- which(!is.na(DcmpRateBIOFv263_1))
K_F <- density(DcmpRateBIOFv263_1[goodInds])
x <- c(K_F$x)
y <- c(K_F$y)
plot(0, 0, col = "white", xlim = c(0.3, 1.0), ylim = c(0, 1.2*max(y)), xlab = "Decay Rate", ylab = "Density", main = "F Pool")
x = c(K_F$x, rev(K_F$x), K_F$x[1])
y = c(rep(0, length(K_F$y)), rev(K_F$y), 0)
y <- constrainDensity(x, y, 0.3, 1.0)
polygon(x, y, col = adjustcolor("purple", 0.25), border = "purple")

x_seq <- seq(0.3, 1.0, by = 0.001)
y_seq <- dnorm(x_seq, 0.66, 0.033)/(1 - pnorm(0.3, 0.66, 0.033))
prior_x <- c(x_seq, rev(x_seq), x_seq[1])
prior_y <- c(y_seq, rep(0, length(y_seq)), y_seq[1])
polygon(prior_x, prior_y, col = adjustcolor("grey", 0.25), border = "grey")
dev.off()



pdf("DecayRate_S.pdf", height = 4, width = 4)
goodInds <- which(!is.na(DcmpRateBIOSv263_1))
K_S <- density(DcmpRateBIOSv263_1[goodInds])
x <- c(K_S$x)
y <- c(K_S$y)
plot(0, 0, col = "white", xlim = c(0.3, 1.0), ylim = c(0, 1.2*max(y)), xlab = "Decay Rate", ylab = "Density", main = "S Pool")
x = c(K_S$x, rev(K_S$x), K_S$x[1])
y = c(rep(0, length(K_S$y)), rev(K_S$y), 0)
y <- constrainDensity(x, y, 0.3, 1.0)
polygon(x, y, col = adjustcolor("purple", 0.25), border = "purple")

x_seq <- seq(0.3, 1.0, by = 0.001)
y_seq <- dnorm(x_seq, 0.66, 0.033)/(1 - pnorm(0.3, 0.66, 0.033))
prior_x <- c(x_seq, rev(x_seq), x_seq[1])
prior_y <- c(y_seq, rep(0, length(y_seq)), y_seq[1])
polygon(prior_x, prior_y, col = adjustcolor("grey", 0.25), border = "grey")
dev.off()



pdf("DecayRate_H.pdf", height = 4, width = 4)
goodInds <- which(!is.na(DcmpRateHUM_1))
K_H <- density(DcmpRateHUM_1[goodInds])
x <- c(K_H$x)
y <- c(K_H$y)
plot(0, 0, col = "white", xlim = c(0.005, 0.05), ylim = c(0, 1.2*max(y)), xlab = "Decay Rate", ylab = "Density", main = "H Pool")
x = c(K_H$x, rev(K_H$x), K_H$x[1])
y = c(rep(0, length(K_H$y)), rev(K_H$y), 0)
y <- constrainDensity(x, y, 0.005, 0.1)
polygon(x, y, col = adjustcolor("purple", 0.25), border = "purple")

x_seq <- seq(0.005, 0.05, by = 0.001)
y_seq <- dnorm(x_seq, 0.02, 0.001)/(1 - pnorm(0.005, 0.02, 0.001))
prior_x <- c(x_seq, rev(x_seq), x_seq[1])
prior_y <- c(y_seq, rep(0, length(y_seq)), y_seq[1])
polygon(prior_x, prior_y, col = adjustcolor("grey", 0.25), border = "grey")
dev.off()





# Fractions of soil carbon shunted between pools

pdf("FracManureDPM.pdf", height = 4, width = 4)
goodInds <- which(!is.na(FracManureDPM_1))
this_density <- density(FracManureDPM_1[goodInds])
x <- c(this_density$x)
y <- c(this_density$y)
plot(0, 0, col = "white", xlim = c(0.4, 0.6), ylim = c(0, 1.2*max(y)), xlab = "FracManureDPM", ylab = "Density", main = "")
x = c(this_density$x, rev(this_density$x), this_density$x[1])
y = c(rep(0, length(this_density$y)), rev(this_density$y), 0)
y <- constrainDensity(x, y, 0.0, 1.0)
polygon(x, y, col = adjustcolor("purple", 0.25), border = "purple")

x_seq <- seq(0.0, 1.0, by = 0.001)
y_seq <- dnorm(x_seq, 0.49, 0.01)
prior_x <- c(x_seq, rev(x_seq), x_seq[1])
prior_y <- c(y_seq, rep(0, length(y_seq)), y_seq[1])
polygon(prior_x, prior_y, col = adjustcolor("grey", 0.25), border = "grey")
dev.off()


pdf("FracManureRPM.pdf", height = 4, width = 4)
goodInds <- which(!is.na(FracManureRPM_1))
this_density <- density(FracManureRPM_1[goodInds])
x <- c(this_density$x)
y <- c(this_density$y)
plot(0, 0, col = "white", xlim = c(0.4, 0.6), ylim = c(0, 1.2*max(y)), xlab = "FracManureRPM", ylab = "Density", main = "")
x = c(this_density$x, rev(this_density$x), this_density$x[1])
y = c(rep(0, length(this_density$y)), rev(this_density$y), 0)
y <- constrainDensity(x, y, 0.0, 1.0)
polygon(x, y, col = adjustcolor("purple", 0.25), border = "purple")

x_seq <- seq(0.0, 1.0, by = 0.001)
y_seq <- dnorm(x_seq, 0.49, 0.01)
prior_x <- c(x_seq, rev(x_seq), x_seq[1])
prior_y <- c(y_seq, rep(0, length(y_seq)), y_seq[1])
polygon(prior_x, prior_y, col = adjustcolor("grey", 0.25), border = "grey")
dev.off()



pdf("FracManureBIOF.pdf", height = 4, width = 4)
goodInds <- which(!is.na(FracManureBIOF_1))
this_density <- density(FracManureBIOF_1[goodInds])
x <- c(this_density$x)
y <- c(this_density$y)
plot(0, 0, col = "white", xlim = c(0.0, 0.1), ylim = c(0, 1.2*max(y)), xlab = "FracManureBIOF", ylab = "Density", main = "")
x = c(this_density$x, rev(this_density$x), this_density$x[1])
y = c(rep(0, length(this_density$y)), rev(this_density$y), 0)
y <- constrainDensity(x, y, 0.0, 1.0)
polygon(x, y, col = adjustcolor("purple", 0.25), border = "purple")

x_seq <- seq(0.0, 1.0, by = 0.001)
y_seq <- dnorm(x_seq, 0.0, 0.01)/0.5
prior_x <- c(x_seq, rev(x_seq), x_seq[1])
prior_y <- c(y_seq, rep(0, length(y_seq)), y_seq[1])
polygon(prior_x, prior_y, col = adjustcolor("grey", 0.25), border = "grey")
dev.off()



pdf("FracManureBIOS.pdf", height = 4, width = 4)
goodInds <- which(!is.na(FracManureBIOS_1))
this_density <- density(FracManureBIOS_1[goodInds])
x <- c(this_density$x)
y <- c(this_density$y)
plot(0, 0, col = "white", xlim = c(0.0, 0.1), ylim = c(0, 1.2*max(y)), xlab = "FracManureBIOS", ylab = "Density", main = "")
x = c(this_density$x, rev(this_density$x), this_density$x[1])
y = c(rep(0, length(this_density$y)), rev(this_density$y), 0)
y <- constrainDensity(x, y, 0.0, 1.0)
polygon(x, y, col = adjustcolor("purple", 0.25), border = "purple")

x_seq <- seq(0.0, 1.0, by = 0.001)
y_seq <- dnorm(x_seq, 0.0, 0.01)/0.5
prior_x <- c(x_seq, rev(x_seq), x_seq[1])
prior_y <- c(y_seq, rep(0, length(y_seq)), y_seq[1])
polygon(prior_x, prior_y, col = adjustcolor("grey", 0.25), border = "grey")
dev.off()




pdf("FracManureHUM.pdf", height = 4, width = 4)
goodInds <- which(!is.na(FracManureHUM_1))
this_density <- density(FracManureHUM_1[goodInds])
x <- c(this_density$x)
y <- c(this_density$y)
plot(0, 0, col = "white", xlim = c(0.0, 0.1), ylim = c(0, 1.2*max(y)), xlab = "FracManureHUM", ylab = "Density", main = "")
x = c(this_density$x, rev(this_density$x), this_density$x[1])
y = c(rep(0, length(this_density$y)), rev(this_density$y), 0)
y <- constrainDensity(x, y, 0.0, 1.0)
polygon(x, y, col = adjustcolor("purple", 0.25), border = "purple")

x_seq <- seq(0.0, 1.0, by = 0.001)
y_seq <- dnorm(x_seq, 0.02, 0.01)/(1 - pnorm(0.0, 0.02, 0.01))
prior_x <- c(x_seq, rev(x_seq), x_seq[1])
prior_y <- c(y_seq, rep(0, length(y_seq)), y_seq[1])
polygon(prior_x, prior_y, col = adjustcolor("grey", 0.25), border = "grey")
dev.off()



pdf("FracSolidXXXToBIOFv263.pdf", height = 4, width = 4)
goodInds <- which(!is.na(FracSolidXXXToBIOFv263_1))
this_density <- density(FracSolidXXXToBIOFv263_1[goodInds])
x <- c(this_density$x)
y <- c(this_density$y)
plot(0, 0, col = "white", xlim = c(0.4, 0.6), ylim = c(0, 1.2*max(y)), xlab = "FracSolidXXXToBIOFv263", ylab = "Density", main = "")
x = c(this_density$x, rev(this_density$x), this_density$x[1])
y = c(rep(0, length(this_density$y)), rev(this_density$y), 0)
y <- constrainDensity(x, y, 0.0, 1.0)
polygon(x, y, col = adjustcolor("purple", 0.25), border = "purple")

x_seq <- seq(0.0, 1.0, by = 0.001)
y_seq <- dnorm(x_seq, 0.46, 0.01)
prior_x <- c(x_seq, rev(x_seq), x_seq[1])
prior_y <- c(y_seq, rep(0, length(y_seq)), y_seq[1])
polygon(prior_x, prior_y, col = adjustcolor("grey", 0.25), border = "grey")
dev.off()



pdf("FracSolidHUMToBIOSv263.pdf", height = 4, width = 4)
goodInds <- which(!is.na(FracSolidHUMToBIOSv263_1))
this_density <- density(FracSolidHUMToBIOSv263_1[goodInds])
x <- c(this_density$x)
y <- c(this_density$y)
plot(0, 0, col = "white", xlim = c(0.4, 0.6), ylim = c(0, 1.2*max(y)), xlab = "FracSolidHUMToBIOSv263", ylab = "Density", main = "")
x = c(this_density$x, rev(this_density$x), this_density$x[1])
y = c(rep(0, length(this_density$y)), rev(this_density$y), 0)
y <- constrainDensity(x, y, 0.0, 1.0)
polygon(x, y, col = adjustcolor("purple", 0.25), border = "purple")

x_seq <- seq(0.0, 1.0, by = 0.001)
y_seq <- dnorm(x_seq, 0.46, 0.01)
prior_x <- c(x_seq, rev(x_seq), x_seq[1])
prior_y <- c(y_seq, rep(0, length(y_seq)), y_seq[1])
polygon(prior_x, prior_y, col = adjustcolor("grey", 0.25), border = "grey")
dev.off()






#other params

pdf("FracClay.pdf", height = 4, width = 4)
goodInds <- which(!is.na(FracClay_1))
this_density <- density(FracClay_1[goodInds])
x <- c(this_density$x)
y <- c(this_density$y)
plot(0, 0, col = "white", xlim = c(0.0, 0.3), ylim = c(0, 1.2*max(y)), xlab = "FracClay", ylab = "Density", main = "")
x = c(this_density$x, rev(this_density$x), this_density$x[1])
y = c(rep(0, length(this_density$y)), rev(this_density$y), 0)
y <- constrainDensity(x, y, 0.0, 1.0)
polygon(x, y, col = adjustcolor("purple", 0.25), border = "purple")

x_seq <- seq(0.0, 1.0, by = 0.001)
y_seq <- dnorm(x_seq, 0.16, 0.02)
prior_x <- c(x_seq, rev(x_seq), x_seq[1])
prior_y <- c(y_seq, rep(0, length(y_seq)), y_seq[1])
polygon(prior_x, prior_y, col = adjustcolor("grey", 0.25), border = "grey")
dev.off()


pdf("RatioDPMToRPM.pdf", height = 4, width = 4)
goodInds <- which(!is.na(RatioDPMToRPM_1))
this_density <- density(RatioDPMToRPM_1[goodInds])
x <- c(this_density$x)
y <- c(this_density$y)
plot(0, 0, col = "white", xlim = c(0.0, 5.0), ylim = c(0, 1.2*max(y)), xlab = "RatioDPMToRPM", ylab = "Density", main = "")
x = c(this_density$x, rev(this_density$x), this_density$x[1])
y = c(rep(0, length(this_density$y)), rev(this_density$y), 0)
y <- constrainDensity(x, y, 0.0, 5.0)
polygon(x, y, col = adjustcolor("purple", 0.25), border = "purple")

x_seq <- seq(0.0, 5.0, by = 0.001)
y_seq <- dnorm(x_seq, 1.44, 0.5)/(pnorm(5, 1.44, 0.5) - pnorm(0, 1.44, 0.5))
prior_x <- c(x_seq, rev(x_seq), x_seq[1])
prior_y <- c(y_seq, rep(0, length(y_seq)), y_seq[1])
polygon(prior_x, prior_y, col = adjustcolor("grey", 0.25), border = "grey")
dev.off()






#Process Error Variances

pdf("sigma2_D.pdf", height = 4, width = 4)
goodInds <- which(!is.na(sigma2_D_1))
this_density <- density(sigma2_D_1[goodInds])
x <- c(this_density$x)
y <- c(this_density$y)
upper = 10E-4
plot(0, 0, col = "white", xlim = c(0.0, upper), ylim = c(0, 1.2*max(y)), xlab = "sigma2_D", ylab = "Density", main = "")
x = c(this_density$x, rev(this_density$x), this_density$x[1])
y = c(rep(0, length(this_density$y)), rev(this_density$y), 0)
y <- constrainDensity(x, y, 0.0, Inf)
polygon(x, y, col = adjustcolor("purple", 0.25), border = "purple")

x_seq <- seq(0.0, upper, by = 0.00001)
y_seq <- dinvgamma(x_seq, shape = 403.42, scale = 1/0.3176)
prior_x <- c(x_seq, rev(x_seq), x_seq[1])
prior_y <- c(y_seq, rep(0, length(y_seq)), y_seq[1])
polygon(prior_x, prior_y, col = adjustcolor("grey", 0.25), border = "grey")
dev.off()



pdf("sigma2_R.pdf", height = 4, width = 4)
goodInds <- which(!is.na(sigma2_R_1))
this_density <- density(sigma2_R_1[goodInds])
x <- c(this_density$x)
y <- c(this_density$y)
upper = 12E-4
plot(0, 0, col = "white", xlim = c(0.0, upper), ylim = c(0, 1.2*max(y)), xlab = "sigma2_R", ylab = "Density", main = "")
x = c(this_density$x, rev(this_density$x), this_density$x[1])
y = c(rep(0, length(this_density$y)), rev(this_density$y), 0)
y <- constrainDensity(x, y, 0.0, Inf)
polygon(x, y, col = adjustcolor("purple", 0.25), border = "purple")

x_seq <- seq(0.0, upper, by = 0.00001)
y_seq <- dinvgamma(x_seq, shape = 403.42, scale = 1/0.3176)
prior_x <- c(x_seq, rev(x_seq), x_seq[1])
prior_y <- c(y_seq, rep(0, length(y_seq)), y_seq[1])
polygon(prior_x, prior_y, col = adjustcolor("grey", 0.25), border = "grey")
dev.off()




pdf("sigma2_F.pdf", height = 4, width = 4)
goodInds <- which(!is.na(sigma2_F_1))
this_density <- density(sigma2_F_1[goodInds])
x <- c(this_density$x)
y <- c(this_density$y)
upper = 10E-4
plot(0, 0, col = "white", xlim = c(0.0, upper), ylim = c(0, 1.2*max(y)), xlab = "sigma2_F", ylab = "Density", main = "")
x = c(this_density$x, rev(this_density$x), this_density$x[1])
y = c(rep(0, length(this_density$y)), rev(this_density$y), 0)
y <- constrainDensity(x, y, 0.0, Inf)
polygon(x, y, col = adjustcolor("purple", 0.25), border = "purple")

x_seq <- seq(0.0, upper, by = 0.00001)
y_seq <- dinvgamma(x_seq, shape = 403.42, scale = 1/0.3176)
prior_x <- c(x_seq, rev(x_seq), x_seq[1])
prior_y <- c(y_seq, rep(0, length(y_seq)), y_seq[1])
polygon(prior_x, prior_y, col = adjustcolor("grey", 0.25), border = "grey")
dev.off()




pdf("sigma2_S.pdf", height = 4, width = 4)
goodInds <- which(!is.na(sigma2_S_1))
this_density <- density(sigma2_S_1[goodInds])
x <- c(this_density$x)
y <- c(this_density$y)
upper = 10E-4
plot(0, 0, col = "white", xlim = c(0.0, upper), ylim = c(0, 1.2*max(y)), xlab = "sigma2_S", ylab = "Density", main = "")
x = c(this_density$x, rev(this_density$x), this_density$x[1])
y = c(rep(0, length(this_density$y)), rev(this_density$y), 0)
y <- constrainDensity(x, y, 0.0, Inf)
polygon(x, y, col = adjustcolor("purple", 0.25), border = "purple")

x_seq <- seq(0.0, upper, by = 0.00001)
y_seq <- dinvgamma(x_seq, shape = 403.42, scale = 1/0.3176)
prior_x <- c(x_seq, rev(x_seq), x_seq[1])
prior_y <- c(y_seq, rep(0, length(y_seq)), y_seq[1])
polygon(prior_x, prior_y, col = adjustcolor("grey", 0.25), border = "grey")
dev.off()



pdf("sigma2_H.pdf", height = 4, width = 4)
goodInds <- which(!is.na(sigma2_H_1))
this_density <- density(sigma2_H_1[goodInds])
x <- c(this_density$x)
y <- c(this_density$y)
upper = 10E-4
plot(0, 0, col = "white", xlim = c(0.0, upper), ylim = c(0, 1.2*max(y)), xlab = "sigma2_H", ylab = "Density", main = "")
x = c(this_density$x, rev(this_density$x), this_density$x[1])
y = c(rep(0, length(this_density$y)), rev(this_density$y), 0)
y <- constrainDensity(x, y, 0.0, Inf)
polygon(x, y, col = adjustcolor("purple", 0.25), border = "purple")

x_seq <- seq(0.0, upper, by = 0.00001)
y_seq <- dinvgamma(x_seq, shape = 403.42, scale = 1/0.3176)
prior_x <- c(x_seq, rev(x_seq), x_seq[1])
prior_y <- c(y_seq, rep(0, length(y_seq)), y_seq[1])
polygon(prior_x, prior_y, col = adjustcolor("grey", 0.25), border = "grey")
dev.off()




#Measurement Error Variances

pdf("sigma2_RPM.pdf", height = 4, width = 4)
goodInds <- which(!is.na(sigma2_RPM_1))
this_density <- density(sigma2_RPM_1[goodInds])
x <- c(this_density$x)
y <- c(this_density$y)
upper = 0.1
plot(0, 0, col = "white", xlim = c(0.0, upper), ylim = c(0, 4*max(y)), xlab = "sigma2_RPM", ylab = "Density", main = "")
x = c(this_density$x, rev(this_density$x), this_density$x[1])
y = c(rep(0, length(this_density$y)), rev(this_density$y), 0)
y <- constrainDensity(x, y, 0.0, Inf)
polygon(x, y, col = adjustcolor("purple", 0.25), border = "purple")

x_seq <- seq(0.0, upper, by = 0.00001)
y_seq <- dinvgamma(x_seq, shape = 10.5, scale =  1/0.03939044)
prior_x <- c(x_seq, rev(x_seq), x_seq[1])
prior_y <- c(y_seq, rep(0, length(y_seq)), y_seq[1])
polygon(prior_x, prior_y, col = adjustcolor("grey", 0.25), border = "grey")
dev.off()



pdf("sigma2_ROC.pdf", height = 4, width = 4)
goodInds <- which(!is.na(sigma2_RPM_1))
this_density <- density(sigma2_RPM_1[goodInds])
x <- c(this_density$x)
y <- c(this_density$y)
upper = 0.1
plot(0, 0, col = "white", xlim = c(0.0, upper), ylim = c(0, 1.2*max(y)), xlab = "sigma2_ROC", ylab = "Density", main = "")
x = c(this_density$x, rev(this_density$x), this_density$x[1])
y = c(rep(0, length(this_density$y)), rev(this_density$y), 0)
y <- constrainDensity(x, y, 0.0, Inf)
polygon(x, y, col = adjustcolor("purple", 0.25), border = "purple")

x_seq <- seq(0.0, upper, by = 0.00001)
y_seq <- dinvgamma(x_seq, shape = 10.5, scale =  1/0.2899684)
prior_x <- c(x_seq, rev(x_seq), x_seq[1])
prior_y <- c(y_seq, rep(0, length(y_seq)), y_seq[1])
polygon(prior_x, prior_y, col = adjustcolor("grey", 0.25), border = "grey")
dev.off()




pdf("sigma2_TOC.pdf", height = 4, width = 4)
goodInds <- which(!is.na(sigma2_TOC_1))
this_density <- density(sigma2_TOC_1[goodInds])
x <- c(this_density$x)
y <- c(this_density$y)
upper = 0.015
plot(0, 0, col = "white", xlim = c(0.0, upper), ylim = c(0, 1.2*max(y)), xlab = "sigma2_TOC", ylab = "Density", main = "")
x = c(this_density$x, rev(this_density$x), this_density$x[1])
y = c(rep(0, length(this_density$y)), rev(this_density$y), 0)
y <- constrainDensity(x, y, 0.0, Inf)
polygon(x, y, col = adjustcolor("purple", 0.25), border = "purple")

x_seq <- seq(0.0, upper, by = 0.00001)
y_seq <- dinvgamma(x_seq, shape = 10.5, scale =   1/0.05292064)
prior_x <- c(x_seq, rev(x_seq), x_seq[1])
prior_y <- c(y_seq, rep(0, length(y_seq)), y_seq[1])
polygon(prior_x, prior_y, col = adjustcolor("grey", 0.25), border = "grey")
dev.off()









# Plot the distributions of the alphas

xlabels <- c("PP", "PF", "Nn0", "Nn1", "Mm0", "Mm1", "Mn0", "Mn1", "In0", "In1", "Im0", "Im1", "Ii0", "Ii1" )

pdf("Treatment_effects.pdf", height = 4, width = 10)
par(mfrow = c(1, 1))
#usedInds <- 469:569
usedInds <- 1:length(treatment_alpha_Fallow_1)
treatment_alpha_Fallow_q <- quantile(exp(treatment_alpha_Fallow_1[usedInds]), probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
treatment_alpha_Pasture_q <- quantile(exp(treatment_alpha_Pasture_1[usedInds]), probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
treatment_alpha_NN0_q <- quantile(exp(treatment_alpha_NN0_1[usedInds]), probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
treatment_alpha_NN1_q <- quantile(exp(treatment_alpha_NN1_1[usedInds]), probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
treatment_alpha_MM0_q <- quantile(exp(treatment_alpha_MM0_1[usedInds]), probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
treatment_alpha_MM1_q <- quantile(exp(treatment_alpha_MM1_1[usedInds]), probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
treatment_alpha_II0_q <- quantile(exp(treatment_alpha_II0_1[usedInds]), probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
treatment_alpha_II1_q <- quantile(exp(treatment_alpha_II1_1[usedInds]), probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
treatment_alpha_IM0_q <- quantile(exp(treatment_alpha_IM0_1[usedInds]), probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
treatment_alpha_IM1_q <- quantile(exp(treatment_alpha_IM1_1[usedInds]), probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
treatment_alpha_IN0_q <- quantile(exp(treatment_alpha_IN0_1[usedInds]), probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
treatment_alpha_IN1_q <- quantile(exp(treatment_alpha_IN1_1[usedInds]), probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
treatment_alpha_MN0_q <- quantile(exp(treatment_alpha_MN0_1[usedInds]), probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
treatment_alpha_MN1_q <- quantile(exp(treatment_alpha_MN1_1[usedInds]), probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)


allQs <- rbind(treatment_alpha_Pasture_q, treatment_alpha_Fallow_q, treatment_alpha_NN0_q, 
treatment_alpha_NN1_q, treatment_alpha_MM0_q, treatment_alpha_MM1_q, treatment_alpha_MN0_q, treatment_alpha_MN1_q, treatment_alpha_IN0_q, 
treatment_alpha_IN1_q, treatment_alpha_IM0_q, treatment_alpha_IM1_q, treatment_alpha_II0_q, 
treatment_alpha_II1_q)



#meds <- allQs[, 3]
#sort_inds <- sort(meds, index.return = TRUE)$ix
sort_inds <- 1:14

plot(0, 0, col = "white", xlim = c(0, 15), ylim = c(0.1, 5), xaxt = "n", xlab = "", ylab = "Treatment Effect", main = "", cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
axis(1, 1:14, xlabels[sort_inds], cex.axis = 1, las = 2, cex.main=1.5, cex.lab=1.5, cex.axis=1.5)

for(i in 1:14)
{
	ix = sort_inds[i]
	lines(x = c(i - 0.25, i + 0.25), y = c(allQs[ix, 3], allQs[ix, 3]), col = adjustcolor("black", 1))
	polygon(x = c(i - 0.25, i + 0.25, i + 0.25, i - 0.25, i - 0.25), y = c(allQs[ix, 2], allQs[ix, 2], allQs[ix, 4], allQs[ix, 4], allQs[ix, 2]), col = adjustcolor("black", 0.25), border = NA)
	polygon(x = c(i - 0.25, i + 0.25, i + 0.25, i - 0.25, i - 0.25), y = c(allQs[ix, 1], allQs[ix, 1], allQs[ix, 5], allQs[ix, 5], allQs[ix, 1]), col = adjustcolor("black", 0.25), border = NA)
}

dev.off()





fieldNames <- c("PP Field 1", "PP Field 2", "PP Field 3",
                "PF Field 1", "PF Field 2", "PF Field 3",
                "MM0 Field 1", "MM0 Field 2", "MM0 Field 3",
                "MM1 Field 1", "MM1 Field 2", "MM1 Field 3",
                "NN0 Field 1", "NN0 Field 2", "NN0 Field 3",
                "NN1 Field 1", "NN1 Field 2", "NN1 Field 3",
                "IM0 Field 1", "IM0 Field 2", "IM0 Field 3",
                "IM1 Field 1", "IM1 Field 2", "IM1 Field 3",
                "IN0 Field 1", "IN0 Field 2", "IN0 Field 3",
                "IN1 Field 1", "IN1 Field 2", "IN1 Field 3",
                "II1 Field 1", "II1 Field 2", "II1 Field 3",
                "II0 Field 1", "II0 Field 2", "II0 Field 3",
                "MN1 Field 1", "MN1 Field 2", "MN1 Field 3",
                "MN0 Field 1", "MN0 Field 2", "MN0 Field 3")



               
treatmentNames <- c("PP",
                "PF",
                "Mm0",
                "Mm1",
                "Nn0",
                "Nn1",
                "Im0",
                "Im1",
                "In0",
                "In1",
                "Ii1",
                "Ii0",
                "Mn1",
                "Mn0")





#Process posterior samples for the latent soil pools in each of the field plots.  To be used for plots.
draws <- draws1

carbonChange <- matrix(NA, 5, 42)
for(f in 1:length(fieldNames))
{
  assign(paste0("D", f), matrix(NA, 5, 113))
  assign(paste0("R", f), matrix(NA, 5, 113))
  assign(paste0("F", f), matrix(NA, 5, 113))
  assign(paste0("S", f), matrix(NA, 5, 113))
  assign(paste0("H", f), matrix(NA, 5, 113))
  assign(paste0("I", f), matrix(NA, 5, 113))
  
  initialC <- draws[, , paste0("D[",f,",", 1, "]")] + draws[, , paste0("R[",f,",", 1, "]")] + draws[, , paste0("F[",f,",", 1, "]")] + draws[, , paste0("S[",f,",", 1, "]")] + draws[, , paste0("H[",f,",", 1, "]")]
  finalC <- draws[, , paste0("D[",f,",", 113, "]")] + draws[, , paste0("R[",f,",", 113, "]")] + draws[, , paste0("F[",f,",", 113, "]")] + draws[, , paste0("S[",f,",", 113, "]")] + draws[, , paste0("H[",f,",", 113, "]")]
  diffC <- initialC - finalC
  carbonChange[, f] <- quantile(diffC, prob = c(0.05, 0.25, 0.5, 0.75, 0.95))
  
  for(t in 1:113)
  {
    quantileD <- quantile(as.numeric(draws[, , paste0("D[",f,",", t, "]")]), prob = 0.05)
    assign(paste0("D", f), `[<-`(eval(as.name(paste0("D", f))), i = 1, j = t, quantileD))
    
    quantileD <- quantile(as.numeric(draws[, , paste0("D[",f,",", t, "]")]), prob = 0.25)
    assign(paste0("D", f), `[<-`(eval(as.name(paste0("D", f))), i = 2, j = t, quantileD))
    
    quantileD <- quantile(as.numeric(draws[, , paste0("D[",f,",", t, "]")]), prob = 0.5)
    assign(paste0("D", f), `[<-`(eval(as.name(paste0("D", f))), i = 3, j = t, quantileD))
    
    quantileD <- quantile(as.numeric(draws[, , paste0("D[",f,",", t, "]")]), prob = 0.75)
    assign(paste0("D", f), `[<-`(eval(as.name(paste0("D", f))), i = 4, j = t, quantileD))
    
    quantileD <- quantile(as.numeric(draws[, , paste0("D[",f,",", t, "]")]), prob = 0.95)
    assign(paste0("D", f), `[<-`(eval(as.name(paste0("D", f))), i = 5, j = t, quantileD))
    
    
    quantileR <- quantile(as.numeric(draws[, , paste0("R[",f,",", t, "]")]), prob = 0.05)
    assign(paste0("R", f), `[<-`(eval(as.name(paste0("R", f))), i = 1, j = t, quantileR))
    
    quantileR <- quantile(as.numeric(draws[, , paste0("R[",f,",", t, "]")]), prob = 0.25)
    assign(paste0("R", f), `[<-`(eval(as.name(paste0("R", f))), i = 2, j = t, quantileR))
    
    quantileR <- quantile(as.numeric(draws[, , paste0("R[",f,",", t, "]")]), prob = 0.5)
    assign(paste0("R", f), `[<-`(eval(as.name(paste0("R", f))), i = 3, j = t, quantileR))
    
    quantileR <- quantile(as.numeric(draws[, , paste0("R[",f,",", t, "]")]), prob = 0.75)
    assign(paste0("R", f), `[<-`(eval(as.name(paste0("R", f))), i = 4, j = t, quantileR))
    
    quantileR <- quantile(as.numeric(draws[, , paste0("R[",f,",", t, "]")]), prob = 0.95)
    assign(paste0("R", f), `[<-`(eval(as.name(paste0("R", f))), i = 5, j = t, quantileR))
    
    
    quantileF <- quantile(as.numeric(draws[, , paste0("F[",f,",", t, "]")]), prob = 0.05)
    assign(paste0("F", f), `[<-`(eval(as.name(paste0("F", f))), i = 1, j = t, quantileF))
    
    quantileF <- quantile(as.numeric(draws[, , paste0("F[",f,",", t, "]")]), prob = 0.25)
    assign(paste0("F", f), `[<-`(eval(as.name(paste0("F", f))), i = 2, j = t, quantileF))
    
    quantileF <- quantile(as.numeric(draws[, , paste0("F[",f,",", t, "]")]), prob = 0.5)
    assign(paste0("F", f), `[<-`(eval(as.name(paste0("F", f))), i = 3, j = t, quantileF))
    
    quantileF <- quantile(as.numeric(draws[, , paste0("F[",f,",", t, "]")]), prob = 0.75)
    assign(paste0("F", f), `[<-`(eval(as.name(paste0("F", f))), i = 4, j = t, quantileF))
    
    quantileF <- quantile(as.numeric(draws[, , paste0("F[",f,",", t, "]")]), prob = 0.95)
    assign(paste0("F", f), `[<-`(eval(as.name(paste0("F", f))), i = 5, j = t, quantileF))
    
    
    quantileS <- quantile(as.numeric(draws[, , paste0("S[",f,",", t, "]")]), prob = 0.05)
    assign(paste0("S", f), `[<-`(eval(as.name(paste0("S", f))), i = 1, j = t, quantileS))
    
    quantileS <- quantile(as.numeric(draws[, , paste0("S[",f,",", t, "]")]), prob = 0.25)
    assign(paste0("S", f), `[<-`(eval(as.name(paste0("S", f))), i = 2, j = t, quantileS))
    
    quantileS <- quantile(as.numeric(draws[, , paste0("S[",f,",", t, "]")]), prob = 0.5)
    assign(paste0("S", f), `[<-`(eval(as.name(paste0("S", f))), i = 3, j = t, quantileS))
    
    quantileS <- quantile(as.numeric(draws[, , paste0("S[",f,",", t, "]")]), prob = 0.75)
    assign(paste0("S", f), `[<-`(eval(as.name(paste0("S", f))), i = 4, j = t, quantileS))
    
    quantileS <- quantile(as.numeric(draws[, , paste0("S[",f,",", t, "]")]), prob = 0.95)
    assign(paste0("S", f), `[<-`(eval(as.name(paste0("S", f))), i = 5, j = t, quantileS))
    
    
    quantileH <- quantile(as.numeric(draws[, , paste0("H[",f,",", t, "]")]), prob = 0.05)
    assign(paste0("H", f), `[<-`(eval(as.name(paste0("H", f))), i = 1, j = t, quantileH))
    
    quantileH <- quantile(as.numeric(draws[, , paste0("H[",f,",", t, "]")]), prob = 0.25)
    assign(paste0("H", f), `[<-`(eval(as.name(paste0("H", f))), i = 2, j = t, quantileH))
    
    quantileH <- quantile(as.numeric(draws[, , paste0("H[",f,",", t, "]")]), prob = 0.5)
    assign(paste0("H", f), `[<-`(eval(as.name(paste0("H", f))), i = 3, j = t, quantileH))
    
    quantileH <- quantile(as.numeric(draws[, , paste0("H[",f,",", t, "]")]), prob = 0.75)
    assign(paste0("H", f), `[<-`(eval(as.name(paste0("H", f))), i = 4, j = t, quantileH))
    
    quantileH <- quantile(as.numeric(draws[, , paste0("H[",f,",", t, "]")]), prob = 0.95)
    assign(paste0("H", f), `[<-`(eval(as.name(paste0("H", f))), i = 5, j = t, quantileH))
    
    
    quantileI <- quantile(as.numeric(draws[, , paste0("IOMCMInit[", f, "]")]), prob = 0.05)
    assign(paste0("I", f), `[<-`(eval(as.name(paste0("I", f))), i = 1, j = t, quantileI))
    
    quantileI <- quantile(as.numeric(draws[, , paste0("IOMCMInit[", f, "]")]), prob = 0.25)
    assign(paste0("I", f), `[<-`(eval(as.name(paste0("I", f))), i = 2, j = t, quantileI))
    
    quantileI <- quantile(as.numeric(draws[, , paste0("IOMCMInit[", f, "]")]), prob = 0.5)
    assign(paste0("I", f), `[<-`(eval(as.name(paste0("I", f))), i = 3, j = t, quantileI))
    
    quantileI <- quantile(as.numeric(draws[, , paste0("IOMCMInit[", f, "]")]), prob = 0.75)
    assign(paste0("I", f), `[<-`(eval(as.name(paste0("I", f))), i = 4, j = t, quantileI))
    
    quantileI <- quantile(as.numeric(draws[, , paste0("IOMCMInit[", f, "]")]), prob = 0.95)
    assign(paste0("I", f), `[<-`(eval(as.name(paste0("I", f))), i = 5, j = t, quantileI))
  }
}

#Rescale to be the carbon change per year
carbonChange <- carbonChange*12/113



#Process posterior samples for the aggregated soil pools acorresponding to each measurement type for each of the field plots.  To be used for plots.

for(f in 1:length(fieldNames))
{
  ROC_temp <- matrix(NA, 5, 113)
  POC_temp <- matrix(NA, 5, 113)
  TOC_temp <- matrix(NA, 5, 113)
  
  ROC_temp_samples <- matrix(NA, 5, 113)
  POC_temp_samples <- matrix(NA, 5, 113)
  TOC_temp_samples <- matrix(NA, 5, 113)
  
  for(t in 1:113)
  {
    thisTOC_samples <- c()
    thisPOC_samples <- c()
    thisROC_samples <- c()
    for(s in 1:100)
    {
      thisTOC_process <- as.numeric(draws[, , paste0("D[",f,",", t, "]")]) + as.numeric(draws[, , paste0("R[",f,",", t, "]")]) + as.numeric(draws[, , paste0("F[",f,",", t, "]")]) + as.numeric(draws[, , paste0("S[",f,",", t, "]")]) + + as.numeric(draws[, , paste0("H[",f,",", t, "]")]) + as.numeric(draws[, , paste0("IOMCMInit[", f, "]")])
      thisPOC_process <- as.numeric(draws[, , paste0("D[",f,",", t, "]")]) + as.numeric(draws[, , paste0("R[",f,",", t, "]")]) + as.numeric(draws[, , paste0("F[",f,",", t, "]")])
      thisROC_process <- as.numeric(draws[, , paste0("IOMCMInit[", f, "]")])
      
      thisTOC_samples <- c(thisTOC_samples, rlnorm(length(thisTOC_process), log(thisTOC_process) - 0.5*sigma2_TOC_1, sqrt(sigma2_TOC_1)))
      thisPOC_samples <- c(thisPOC_samples, rlnorm(length(thisPOC_process), log(thisPOC_process) - 0.5*sigma2_RPM_1, sqrt(sigma2_RPM_1)))
      thisROC_samples <- c(thisROC_samples, rlnorm(length(thisROC_process), log(thisROC_process) - 0.5*sigma2_ROC_1, sqrt(sigma2_ROC_1)))
    }
    
    ROC_temp[, t] <- quantile(thisROC_process, prob = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
    POC_temp[, t] <- quantile(thisPOC_process, prob = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
    TOC_temp[, t] <- quantile(thisTOC_process, prob = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
    
    ROC_temp_samples[, t] <- quantile(thisROC_samples, prob = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
    POC_temp_samples[, t] <- quantile(thisPOC_samples, prob = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
    TOC_temp_samples[, t] <- quantile(thisTOC_samples, prob = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
  }
  
  assign(paste0("ROC", f), ROC_temp)
  assign(paste0("POC", f), POC_temp)
  assign(paste0("TOC", f), TOC_temp)
  
  assign(paste0("ROC_plusNoise", f), ROC_temp_samples)
  assign(paste0("POC_plusNoise", f), POC_temp_samples)
  assign(paste0("TOC_plusNoise", f), TOC_temp_samples)
}






#Create plots for each field-plot for each of the latent soil carbon pools.
for(f in 1:length(fieldNames))
{
  tiff(file = paste0(fieldNames[f], ".tif"), width = 8, height = 16, units = "in", res = 200)
  
  thisD <- get(paste0("D", f))
  thisR <- get(paste0("R", f))
  thisF <- get(paste0("F", f))
  thisS <- get(paste0("S", f))
  thisH <- get(paste0("H", f))
  thisI <- get(paste0("I", f))
  
  par(mfrow = c(6, 1), mar = c(5.1, 5.1, 4.1, 2.1), oma = c(0,0,0,0))
  plot(x = c(0, 113), y = c(0, 0), col = adjustcolor("white", 0), xlim = c(0, 113), ylim = c(0, max(thisD)), xlab = "Months", ylab = "Carbon (t/ha)", main = "D Pool", cex.main = 2.2, cex.lab = 2.2, cex.axis = 2.2)
  polygon(x = c(1:113, rev(1:113), 1), y = c(thisD[1, ], rev(thisD[5, ]), thisD[1, 1]), col = adjustcolor("blue", 0.25), border = NA)
  polygon(x = c(1:113, rev(1:113), 1), y = c(thisD[2, ], rev(thisD[4, ]), thisD[2, 1]), col = adjustcolor("blue", 0.25), border = NA)
  lines(x = 1:113, y = thisD[3, ], col = "blue")
  
  plot(x = c(0, 113), y = c(0, 0), col = adjustcolor("white", 0), xlim = c(0, 113), ylim = c(0, max(thisR)), xlab = "Months", ylab = "Carbon (t/ha)", main = "R Pool", cex.main = 2.2, cex.lab = 2.2, cex.axis = 2.2)
  polygon(x = c(1:113, rev(1:113), 1), y = c(thisR[1, ], rev(thisR[5, ]), thisR[1, 1]), col = adjustcolor("blue", 0.25), border = NA)
  polygon(x = c(1:113, rev(1:113), 1), y = c(thisR[2, ], rev(thisR[4, ]), thisR[2, 1]), col = adjustcolor("blue", 0.25), border = NA)
  lines(x = 1:113, y = thisR[3, ], col = "blue")
  
  plot(x = c(0, 113), y = c(0, 0), col = adjustcolor("white", 0), xlim = c(0, 113), ylim = c(0, max(thisF)), xlab = "Months", ylab = "Carbon (t/ha)", main = "F Pool", cex.main = 2.2, cex.lab = 2.2, cex.axis = 2.2)
  polygon(x = c(1:113, rev(1:113), 1), y = c(thisF[1, ], rev(thisF[5, ]), thisF[1, 1]), col = adjustcolor("blue", 0.25), border = NA)
  polygon(x = c(1:113, rev(1:113), 1), y = c(thisF[2, ], rev(thisF[4, ]), thisF[2, 1]), col = adjustcolor("blue", 0.25), border = NA)
  lines(x = 1:113, y = thisF[3, ], col = "blue")
  
  plot(x = c(0, 113), y = c(0, 0), col = adjustcolor("white", 0), xlim = c(0, 113), ylim = c(0, max(thisS)), xlab = "Months", ylab = "Carbon (t/ha)", main = "S Pool", cex.main = 2.2, cex.lab = 2.2, cex.axis = 2.2)
  polygon(x = c(1:113, rev(1:113), 1), y = c(thisS[1, ], rev(thisS[5, ]), thisS[1, 1]), col = adjustcolor("blue", 0.25), border = NA)
  polygon(x = c(1:113, rev(1:113), 1), y = c(thisS[2, ], rev(thisS[4, ]), thisS[2, 1]), col = adjustcolor("blue", 0.25), border = NA)
  lines(x = 1:113, y = thisS[3, ], col = "blue")
  
  plot(x = c(0, 113), y = c(0, 0), col = adjustcolor("white", 0), xlim = c(0, 113), ylim = c(0, max(thisH)), xlab = "Months", ylab = "Carbon (t/ha)", main = "H Pool", cex.main = 2.2, cex.lab = 2.2, cex.axis = 2.2)
  polygon(x = c(1:113, rev(1:113), 1), y = c(thisH[1, ], rev(thisH[5, ]), thisH[1, 1]), col = adjustcolor("blue", 0.25), border = NA)
  polygon(x = c(1:113, rev(1:113), 1), y = c(thisH[2, ], rev(thisH[4, ]), thisH[2, 1]), col = adjustcolor("blue", 0.25), border = NA)
  lines(x = 1:113, y = thisH[3, ], col = "blue")
  
  plot(x = c(0, 113), y = c(0, 0), col = adjustcolor("white", 0), xlim = c(0, 113), ylim = c(0, max(thisI)), xlab = "Months", ylab = "Carbon (t/ha)", main = "I Pool", cex.main = 2.2, cex.lab = 2.2, cex.axis = 2.2)
  polygon(x = c(1:113, rev(1:113), 1), y = c(thisI[1, ], rev(thisI[5, ]), thisI[1, 1]), col = adjustcolor("blue", 0.25), border = NA)
  polygon(x = c(1:113, rev(1:113), 1), y = c(thisI[2, ], rev(thisI[4, ]), thisI[2, 1]), col = adjustcolor("blue", 0.25), border = NA)
  lines(x = 1:113, y = thisI[3, ], col = "blue")
  
  dev.off()
}


#Create plots for each field-plot for each of the latent soil carbon pools.
pdf(file = "Latent_Pools.pdf", width = 10, height = 20)
for(f in 1:length(fieldNames))
{
  thisD <- get(paste0("D", f))
  thisR <- get(paste0("R", f))
  thisF <- get(paste0("F", f))
  thisS <- get(paste0("S", f))
  thisH <- get(paste0("H", f))
  thisI <- get(paste0("I", f))
  
  par(mfrow = c(6, 1), mar = c(5.1, 5.1, 4.1, 2.1), oma = c(0,0,0,0))
  plot(x = c(0, 113), y = c(0, 0), col = adjustcolor("white", 0), xlim = c(0, 113), ylim = c(0, max(thisD)), xlab = "Months", ylab = "Carbon (t/ha)", main = paste0("D Pool (", fieldNames[f], ")"), cex.main = 2.2, cex.lab = 2.2, cex.axis = 2.2)
  polygon(x = c(1:113, rev(1:113), 1), y = c(thisD[1, ], rev(thisD[5, ]), thisD[1, 1]), col = adjustcolor("blue", 0.25), border = NA)
  polygon(x = c(1:113, rev(1:113), 1), y = c(thisD[2, ], rev(thisD[4, ]), thisD[2, 1]), col = adjustcolor("blue", 0.25), border = NA)
  lines(x = 1:113, y = thisD[3, ], col = "blue")
  
  plot(x = c(0, 113), y = c(0, 0), col = adjustcolor("white", 0), xlim = c(0, 113), ylim = c(0, max(thisR)), xlab = "Months", ylab = "Carbon (t/ha)", main = paste0("R Pool: ", fieldNames[f]), cex.main = 2.2, cex.lab = 2.2, cex.axis = 2.2)
  polygon(x = c(1:113, rev(1:113), 1), y = c(thisR[1, ], rev(thisR[5, ]), thisR[1, 1]), col = adjustcolor("blue", 0.25), border = NA)
  polygon(x = c(1:113, rev(1:113), 1), y = c(thisR[2, ], rev(thisR[4, ]), thisR[2, 1]), col = adjustcolor("blue", 0.25), border = NA)
  lines(x = 1:113, y = thisR[3, ], col = "blue")
  
  plot(x = c(0, 113), y = c(0, 0), col = adjustcolor("white", 0), xlim = c(0, 113), ylim = c(0, max(thisF)), xlab = "Months", ylab = "Carbon (t/ha)", main = paste0("F Pool: ", fieldNames[f]), cex.main = 2.2, cex.lab = 2.2, cex.axis = 2.2)
  polygon(x = c(1:113, rev(1:113), 1), y = c(thisF[1, ], rev(thisF[5, ]), thisF[1, 1]), col = adjustcolor("blue", 0.25), border = NA)
  polygon(x = c(1:113, rev(1:113), 1), y = c(thisF[2, ], rev(thisF[4, ]), thisF[2, 1]), col = adjustcolor("blue", 0.25), border = NA)
  lines(x = 1:113, y = thisF[3, ], col = "blue")
  
  plot(x = c(0, 113), y = c(0, 0), col = adjustcolor("white", 0), xlim = c(0, 113), ylim = c(0, max(thisS)), xlab = "Months", ylab = "Carbon (t/ha)", main = paste0("S Pool: ", fieldNames[f]), cex.main = 2.2, cex.lab = 2.2, cex.axis = 2.2)
  polygon(x = c(1:113, rev(1:113), 1), y = c(thisS[1, ], rev(thisS[5, ]), thisS[1, 1]), col = adjustcolor("blue", 0.25), border = NA)
  polygon(x = c(1:113, rev(1:113), 1), y = c(thisS[2, ], rev(thisS[4, ]), thisS[2, 1]), col = adjustcolor("blue", 0.25), border = NA)
  lines(x = 1:113, y = thisS[3, ], col = "blue")
  
  plot(x = c(0, 113), y = c(0, 0), col = adjustcolor("white", 0), xlim = c(0, 113), ylim = c(0, max(thisH)), xlab = "Months", ylab = "Carbon (t/ha)", main = paste0("H Pool: ", fieldNames[f]), cex.main = 2.2, cex.lab = 2.2, cex.axis = 2.2)
  polygon(x = c(1:113, rev(1:113), 1), y = c(thisH[1, ], rev(thisH[5, ]), thisH[1, 1]), col = adjustcolor("blue", 0.25), border = NA)
  polygon(x = c(1:113, rev(1:113), 1), y = c(thisH[2, ], rev(thisH[4, ]), thisH[2, 1]), col = adjustcolor("blue", 0.25), border = NA)
  lines(x = 1:113, y = thisH[3, ], col = "blue")
  
  plot(x = c(0, 113), y = c(0, 0), col = adjustcolor("white", 0), xlim = c(0, 113), ylim = c(0, max(thisI)), xlab = "Months", ylab = "Carbon (t/ha)", main = paste0("I Pool: ", fieldNames[f]), cex.main = 2.2, cex.lab = 2.2, cex.axis = 2.2)
  polygon(x = c(1:113, rev(1:113), 1), y = c(thisI[1, ], rev(thisI[5, ]), thisI[1, 1]), col = adjustcolor("blue", 0.25), border = NA)
  polygon(x = c(1:113, rev(1:113), 1), y = c(thisI[2, ], rev(thisI[4, ]), thisI[2, 1]), col = adjustcolor("blue", 0.25), border = NA)
  lines(x = 1:113, y = thisI[3, ], col = "blue")
  
}
dev.off()




draws <- draws1
#Calculate the carbon change for each field-plot
carbonChange_treatments <- matrix(NA, 5, 14)
carbonProbs <- rep(NA, 14)
for(f in 1:14)
{
  assign(paste0("D", f), matrix(NA, 5, 113))
  assign(paste0("R", f), matrix(NA, 5, 113))
  assign(paste0("F", f), matrix(NA, 5, 113))
  assign(paste0("S", f), matrix(NA, 5, 113))
  assign(paste0("H", f), matrix(NA, 5, 113))
  assign(paste0("I", f), matrix(NA, 5, 113))
  
  f_ind <- (f - 1)*3 + 1
  initialC1 <- draws[, , paste0("D[",f_ind,",", 1, "]")] + draws[, , paste0("R[",f_ind,",", 1, "]")] + draws[, , paste0("F[",f_ind,",", 1, "]")] + draws[, , paste0("S[",f_ind,",", 1, "]")] + draws[, , paste0("H[",f_ind,",", 1, "]")]
  finalC1 <- draws[, , paste0("D[",f_ind,",", 113, "]")] + draws[, , paste0("R[",f_ind,",", 113, "]")] + draws[, , paste0("F[",f_ind,",", 113, "]")] + draws[, , paste0("S[",f_ind,",", 113, "]")] + draws[, , paste0("H[",f_ind,",", 113, "]")]
  
  f_ind <- (f - 1)*3 + 2
  initialC2 <- draws[, , paste0("D[",f_ind,",", 1, "]")] + draws[, , paste0("R[",f_ind,",", 1, "]")] + draws[, , paste0("F[",f_ind,",", 1, "]")] + draws[, , paste0("S[",f_ind,",", 1, "]")] + draws[, , paste0("H[",f_ind,",", 1, "]")]
  finalC2 <- draws[, , paste0("D[",f_ind,",", 113, "]")] + draws[, , paste0("R[",f_ind,",", 113, "]")] + draws[, , paste0("F[",f_ind,",", 113, "]")] + draws[, , paste0("S[",f_ind,",", 113, "]")] + draws[, , paste0("H[",f_ind,",", 113, "]")]

  f_ind <- (f - 1)*3 + 3
  initialC3 <- draws[, , paste0("D[",f_ind,",", 1, "]")] + draws[, , paste0("R[",f_ind,",", 1, "]")] + draws[, , paste0("F[",f_ind,",", 1, "]")] + draws[, , paste0("S[",f_ind,",", 1, "]")] + draws[, , paste0("H[",f_ind,",", 1, "]")]
  finalC3 <- draws[, , paste0("D[",f_ind,",", 113, "]")] + draws[, , paste0("R[",f_ind,",", 113, "]")] + draws[, , paste0("F[",f_ind,",", 113, "]")] + draws[, , paste0("S[",f_ind,",", 113, "]")] + draws[, , paste0("H[",f_ind,",", 113, "]")]

  #initialC <- c(initialC1, initialC2, initialC3)
  initialC <- (initialC1 + initialC2 + initialC3)/3
  finalC <- (finalC1 + finalC2 + finalC3)/3
  diffC <- initialC - finalC
  carbonChange_treatments[, f] <- quantile(diffC, prob = c(0.05, 0.25, 0.5, 0.75, 0.95))
  carbonProbs[f] <- round(length(which(diffC > 0))/length(diffC), 2)
}
#Rescale to be the carbon change per year
carbonChange_treatments <- carbonChange_treatments*12/113



#Key to numeric values used in the Stan model for treatments
#1"Pasture",
#2"Fallow",
#3"MM0",
#4"MM1",
#5"NN0",
#6"NN1",
#7"IM0",
#8"IM1",
#9"IN0",
#10"IN1",
#11"II1",
#12"II0",
#13"MN1",
#14"MN0")


#Plot the carbon change by treatments
pdf("carbon_change_byTreatment.pdf", width = 12)
sortInd <- c(1, 2, 5, 6, 3, 4, 14, 13, 9, 10, 7, 8, 12, 11)
fieldCols <- rep("black", 14)
par(mar = c(7, 5, 1, 1))
plot(1:14, rep(0, 14), col = "white", ylim = c(-1.5, 3), xaxt = "n", ylab = "Carbon Flux (t/ha/y)", xlab = "", cex.main = 2, cex.lab = 2, cex.axis = 2)
axis(1, 1:14, treatmentNames[sortInd], cex.axis = 1, las = 2, , cex.main = 2, cex.lab = 2, cex.axis = 2)
for(i in 1:14)
{
  x = c(i - 0.25, i + 0.25, i + 0.25, i - 0.25)
  y1 = c(carbonChange_treatments[1, sortInd[i]], carbonChange_treatments[1, sortInd[i]], carbonChange_treatments[5, sortInd[i]], carbonChange_treatments[5, sortInd[i]])
  y2 = c(carbonChange_treatments[2, sortInd[i]], carbonChange_treatments[2, sortInd[i]], carbonChange_treatments[4, sortInd[i]], carbonChange_treatments[4, sortInd[i]])
  polygon(x = x, y = y1, col = adjustcolor(fieldCols[sortInd[i]], 0.25), border = NA)
  polygon(x = x, y = y2, col = adjustcolor(fieldCols[sortInd[i]], 0.25), border = NA)
  lines(c(i - 0.25, i + 0.25), c(carbonChange_treatments[3, sortInd[i]], carbonChange_treatments[3, sortInd[i]]), col = fieldCols[sortInd[i]])
  text(1 - carbonProbs[sortInd[i]], x = i, y = carbonChange_treatments[5, sortInd[i]] + 0.1, cex = 1.25)
}
lines(c(0, 15), c(0, 0), col = "red", lty = 2)

dev.off()


#Plot the carbon change by cover crop
pdf("carbon_change_byCoverCrop.pdf", width = 6)
fieldCols <- rep("black", 2)
par(mar = c(7, 5, 1, 1))
plot(1:2, rep(0, 2), col = "white", ylim = c(-0.5, 0.5), xaxt = "n", ylab = "Mean Carbon Flux (t/ha/y)", xlab = "", xlim = c(0, 3))
axis(1, 1:2, c("Cover Crop", "No Cover Crop"), cex.axis = 1, las = 1)
for(i in 1:2)
{
  x = c(i - 0.25, i + 0.25, i + 0.25, i - 0.25)
  y1 = c(carbonChange_cover[1, i], carbonChange_cover[1, i], carbonChange_cover[5, i], carbonChange_cover[5, i])
  y2 = c(carbonChange_cover[2, i], carbonChange_cover[2, i], carbonChange_cover[4, i], carbonChange_cover[4, i])
  polygon(x = x, y = y1, col = adjustcolor(fieldCols[i], 0.25), border = NA)
  polygon(x = x, y = y2, col = adjustcolor(fieldCols[i], 0.25), border = NA)
  lines(c(i - 0.25, i + 0.25), c(carbonChange_cover[3, i], carbonChange_cover[3, i]), col = fieldCols[i])
  text(carbonProbs_01[i], x = i, y = carbonChange_cover[5, i] + 0.05)
}
lines(c(0, 15), c(0, 0), col = "red", lty = 2)
dev.off()





#Import the data for the Millennium Tilalge Trial (MTT).
source("data_cmdStan_truncated_allTreatments")

#Plot the measurements with the corresponding combination of latent pools for each field-plot in the MTT.
for(f in 1:length(fieldNames))
{
  thisPOC_plusNoise <- get(paste0("POC_plusNoise", f))
  thisTOC_plusNoise <- get(paste0("TOC_plusNoise", f))
  thisROC_plusNoise <- get(paste0("ROC_plusNoise", f))
  
  thisPOC <- get(paste0("POC", f))
  thisTOC <- get(paste0("TOC", f))
  thisROC <- get(paste0("ROC", f))
  
  tiff(file = paste0(fieldNames[f], "_obs.tif"), width = 8, height = 16, units = "in", res = 200)
  
  par(mfrow = c(3, 1), mar = c(5.1, 5.1, 4.1, 2.1), oma = c(0,0,0,0))
  
  ylim1 <- c(0.9*min(thisPOC_plusNoise), 1.1*max(thisPOC_plusNoise))
  plot(x = c(0, 113), y = c(0, 0), col = adjustcolor("white", 0), xlim = c(0, 113), ylim = ylim1, xlab = "Months", ylab = "Carbon (t/ha)", main = "POC", cex.main = 2, cex.lab = 2, cex.axis = 2)
  polygon(x = c(1:113, rev(1:113), 1), y = c(thisPOC_plusNoise[1, ], rev(thisPOC_plusNoise[5, ]), thisPOC_plusNoise[1, 1]), col = adjustcolor("grey", 0.25), border = NA)
  polygon(x = c(1:113, rev(1:113), 1), y = c(thisPOC_plusNoise[2, ], rev(thisPOC_plusNoise[4, ]), thisPOC_plusNoise[2, 1]), col = adjustcolor("grey", 0.25), border = NA)
  polygon(x = c(1:113, rev(1:113), 1), y = c(thisPOC[1, ], rev(thisPOC[5, ]), thisPOC[1, 1]), col = adjustcolor("purple", 0.25), border = NA)
  polygon(x = c(1:113, rev(1:113), 1), y = c(thisPOC[2, ], rev(thisPOC[4, ]), thisPOC[2, 1]), col = adjustcolor("purple", 0.25), border = NA)
  lines(x = 1:113, y = thisPOC[3, ], col = adjustcolor("purple", 0.5))
  #lines(x = 1:113, y = thisPOC_plusNoise[3, ], col = adjustcolor("grey", 0.5))
  lines(x = 1:113, y = thisPOC_plusNoise[1, ], lty = "dashed")
  lines(x = 1:113, y = thisPOC_plusNoise[5, ], lty = "dashed")
  ind <- which(mPOC[f, ] != -9999)
  points(x = ind, y = mPOC[f, ind], cex = 2, pch = 20)
  
  ylim2 <- c(0.9*min(thisROC_plusNoise), 1.1*max(thisROC_plusNoise))
  plot(x = c(0, 113), y = c(0, 0), col = adjustcolor("white", 0), ylim = ylim2, xlim = c(0, 113), xlab = "Months", ylab = "Carbon (t/ha)", main = "ROC", cex.main = 2, cex.lab = 2, cex.axis = 2)
  polygon(x = c(1:113, rev(1:113), 1), y = c(thisROC_plusNoise[1, ], rev(thisROC_plusNoise[5, ]), thisROC_plusNoise[1, 1]), col = adjustcolor("grey", 0.25), border = NA)
  polygon(x = c(1:113, rev(1:113), 1), y = c(thisROC_plusNoise[2, ], rev(thisROC_plusNoise[4, ]), thisROC_plusNoise[2, 1]), col = adjustcolor("grey", 0.25), border = NA)
  polygon(x = c(1:113, rev(1:113), 1), y = c(thisROC[1, ], rev(thisROC[5, ]), thisROC[1, 1]), col = adjustcolor("purple", 0.25), border = NA)
  polygon(x = c(1:113, rev(1:113), 1), y = c(thisROC[2, ], rev(thisROC[4, ]), thisROC[2, 1]), col = adjustcolor("purple", 0.25), border = NA)
  lines(x = 1:113, y = thisROC[3, ], col = adjustcolor("purple", 0.5))
  #lines(x = 1:113, y = thisROC_plusNoise[3, ], col = adjustcolor("grey", 0.5))
  lines(x = 1:113, y = thisROC_plusNoise[1, ], lty = "dashed")
  lines(x = 1:113, y = thisROC_plusNoise[5, ], lty = "dashed")
  ind <- which(mROC[f, ] != -9999)
  points(x = ind, y = mROC[f, ind], cex = 2, pch = 20)
  
  ylim3 <- c(0.9*min(thisTOC_plusNoise), 1.1*max(thisTOC_plusNoise))
  plot(x = c(0, 113), y = c(0, 0), col = adjustcolor("white", 0), xlim = c(0, 113), ylim = ylim3, xlab = "Months", ylab = "Carbon (t/ha)", main = "TOC", cex.main = 2, cex.lab = 2, cex.axis = 2)
  polygon(x = c(1:113, rev(1:113), 1), y = c(thisTOC_plusNoise[1, ], rev(thisTOC_plusNoise[5, ]), thisTOC_plusNoise[1, 1]), col = adjustcolor("grey", 0.25), border = NA)
  polygon(x = c(1:113, rev(1:113), 1), y = c(thisTOC_plusNoise[2, ], rev(thisTOC_plusNoise[4, ]), thisTOC_plusNoise[2, 1]), col = adjustcolor("grey", 0.25), border = NA)
  polygon(x = c(1:113, rev(1:113), 1), y = c(thisTOC[1, ], rev(thisTOC[5, ]), thisTOC[1, 1]), col = adjustcolor("purple", 0.25), border = NA)
  polygon(x = c(1:113, rev(1:113), 1), y = c(thisTOC[2, ], rev(thisTOC[4, ]), thisTOC[2, 1]), col = adjustcolor("purple", 0.25), border = NA)
  lines(x = 1:113, y = thisTOC[3, ], col = adjustcolor("purple", 0.5))
  #lines(x = 1:113, y = thisTOC_plusNoise[3, ], col = adjustcolor("grey", 0.5))
  lines(x = 1:113, y = thisTOC_plusNoise[1, ], lty = "dashed")
  lines(x = 1:113, y = thisTOC_plusNoise[5, ], lty = "dashed")
  ind <- which(mTOC[f, ] != -9999)
  points(x = ind, y = mTOC[f, ind], cex = 2, pch = 20)
  dev.off()
}



