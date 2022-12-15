# call packages -----------------------------------------------------------
packages <- c("openxlsx",
              "ggplot2",
              "tidymv",
              "car",
              "ggfortify",
              "gridExtra",
              "grid",
              "ggpubr",
              "openair",
              "cvTools",
              "extrafont",
              "remotes",
              "R2jags",
              "dplyr",
              "broom.mixed",
              "boot",
              "ggmcmc",
              "scales",
              "postpack",
              "MCMCvis")

if (!require(install.load)) {
  install.packages("install.load")
}

install.load::install_load(packages)

# reconcile fonts ---------------------------------------------------------
remotes::install_version("Rttf2pt1", version = "1.3.8")
font_import(prompt = FALSE, pattern = "calibri")
fonts()
loadfonts(device = "win")
windowsFonts()

# data steps --------------------------------------------------------------
## read-in raw data
cbpSARequiv.dat <- read.xlsx("Input\\SARvLWGsp-su_esc.xlsx",
                             sheet = 1,
                             colNames = TRUE)

## create reduced data set
cbpSARequiv.datRed <- subset(cbpSARequiv.dat, select = c(src,sar,esc))
cbpSARequiv.outwb <- createWorkbook()
addWorksheet(cbpSARequiv.outwb, "lawry")
addWorksheet(cbpSARequiv.outwb, "Raymond")

## define CBP benchmarks
lowCBP <- 43000
medCBP <- 137000
hiCBP <- 235000

# lawry analysis ----------------------------------------------------------

## define initial and other values
initLaw.mod <- lm(data=subset(cbpSARequiv.dat,src=="lawry"),sar~esc)
initLaw.alpha <- as.numeric(initLaw.mod$coef[1])
initLaw.alphaSE <- coef(summary(initLaw.mod))[1,2]
initLaw.beta <- as.numeric(initLaw.mod$coef[2])
initLaw.betaSE <- coef(summary(initLaw.mod))[2,2]
df.law <- summary(initLaw.mod)$df[2]

## Fit model
### formulate data input
cbpSARequiv.lawDat <- list(esc.law = subset(cbpSARequiv.dat,src=="lawry")[,3],
                           sar.law = subset(cbpSARequiv.dat,src=="lawry")[,4],
                           n.law = nrow(subset(cbpSARequiv.dat,src=="lawry")))

### define initial values
inits.law <- function()
{
  list(alpha.law = runif(1,initLaw.alpha-(5*initLaw.alphaSE),initLaw.alpha+(5*initLaw.alphaSE)),
       beta.law = runif(1,initLaw.beta-(5*initLaw.betaSE),initLaw.beta+(5*initLaw.betaSE)))
}

### specify model
cat('

model {
  # Likelihood
  for (i in 1:n.law){
    sar.law[i] ~ dt(mu.law[i],tau.law,nu.law)
    mu.law[i] <- alpha.law + beta.law * esc.law[i]

    # # res.law[i] <- sar.law[i] - mu.law[i]
    # # sar.new[i] ~ dt(mu.law[i],tau.law,nu.law)
    # # res.sar.new[i] <- sar.new[i] - mu.law[i]
    #
    # posterior predictive check
    ## predicted values
    law.pred[i]<-mu.law[i]
    ## residuals for observed data
    law.resid[i]<-sar.law[i]-law.pred[i]
    ## discrepancy/squared residuals for observed data
    law.SqResid[i]<-pow(law.resid[i],2)

    # generate replicate data and compute fit stats
    sar.new[i]~dt(mu.law[i],tau.law,nu.law)
    # discrepancy/squared residuals for new/ideal data
    law.SqResid.new[i]<-pow(sar.new[i]-law.pred[i],2)
  }

  # Priors
  alpha.law ~ dnorm(0 , 0.001)
  beta.law ~ dnorm(0 , 0.001)
  nuMinusOne.law ~ dexp(1/29)
  nu.law <- nuMinusOne.law+1
  tau.law <- 1 / (sigma.law * sigma.law)
  sigma.law~dunif(0,100)

  # predict sar at different levels of CBP goals
  ## define CBP levels
  lowCBP <- 43000
  medCBP <- 137000
  hiCBP <- 235000

  ## generate predictions
  muLow <- alpha.law + beta.law*lowCBP
  muMed <- alpha.law + beta.law*medCBP
  muHi <- alpha.law + beta.law*hiCBP
  sarLow ~ dt(muLow,tau.law,nu.law)
  sarMed ~ dt(muMed,tau.law,nu.law)
  sarHi ~ dt(muHi,tau.law,nu.law)

  # derived values
  ## sum of squared residuals for actual dataset
  law.sumsqO<-sum(law.SqResid[])
  ## sum of squared residuals for new dataset
  law.sumsqN<-sum(law.SqResid.new[])
  ## test whether new data are more extreme
  law.test<-step(law.sumsqN-law.sumsqO)
  ## bayesian pvalue
  law.bpval<-mean(law.test)


}', file={lmLaw.jags <- tempfile()})

### define parameters to monitor
params.law <- c("alpha.law",
                "beta.law",
                "sigma.law",
                "nu.law",
                "sarLow",
                "sarMed",
                "sarHi",
                "law.sumsqO",
                "law.sumsqN",
                "law.bpval",
                "sar.new")

### call jags
fit.law <- jags(data = cbpSARequiv.lawDat,
                inits = inits.law,
                parameters.to.save = params.law,
                model.file = lmLaw.jags,
                n.chains = 3,
                n.iter = 250000,
                n.burnin = 25000,
                n.thin = 10,
                DIC = F)

## review output
fit.law

## review diagnostics
### format output
law.mcmc <- as.mcmc(fit.law)
law.ggs.dat <- ggs(law.mcmc)

law.ppc.dat <- data.frame(SumSqO=as.data.frame(subset(law.ggs.dat,Parameter=="law.sumsqO",select=c(value))),
                          SumSqN=as.data.frame(subset(law.ggs.dat,Parameter=="law.sumsqN",select=c(value))))

colnames(law.ppc.dat) <- c("law.sumsqO",
                           "law.sumsqN")
head(law.ppc.dat)

max(law.ppc.dat$law.sumsqN)

### trace plots
#### alpha
tracePlot.alpha.law <- ggs_traceplot(law.ggs.dat, family = "alpha.law")+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title.y = element_text(face = "bold", size = 16,vjust = 1,margin = margin(t = 0, r = 50, b = 0, l = 0),family = "Calibri"),
        axis.title.x = element_text(face = "bold", size = 16,vjust = -1,margin = margin(t = 10, r = 0, b = 0, l = 0),family = "Calibri"),
        axis.text.x = element_text(face = "bold",size = 14,color="black", vjust=0.5,family = "Calibri"),
        axis.text.y = element_text(face = "bold",size = 14,color="black",family = "Calibri"),
        strip.background = element_blank(),
        legend.position = "none",
        strip.text.x = element_blank(),
        legend.title = element_text(face = "bold",size = 12,color="black",family = "Calibri"),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        legend.text=element_text(face = "bold",size = 12,color="black",family = "Calibri"),
        axis.ticks.length = unit(0.15, "cm"))+
  labs(title = "alpha (Lawry [2020] & McCann [2022])",y = "Value", x = "Iteration") +
  annotate(geom = "text",
           x = post_dim(law.mcmc,"saved")/post_dim(law.mcmc,"chains"),
           y = max(MCMCchains(law.mcmc, params = 'alpha.law')),
           label = paste("hat(R)","~`=`~",round(as.numeric(fit.law$BUGSoutput$summary[1,8]),3),sep=""),
           hjust = 1.0,vjust=0.0,size = 5.1,family = "Calibri",parse = TRUE)+ # Rhat has to be Changed based on jags output (extraction rounds to nearest interger)
  theme(plot.title = element_text(hjust = 0.5,size = 16,face = "bold",family = "Calibri"))

#### beta
tracePlot.beta.law <- ggs_traceplot(law.ggs.dat, family = "beta.law")+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title.y = element_text(face = "bold", size = 16,vjust = 1,margin = margin(t = 0, r = 50, b = 0, l = 0),family = "Calibri"),
        axis.title.x = element_text(face = "bold", size = 16,vjust = -1,margin = margin(t = 10, r = 0, b = 0, l = 0),family = "Calibri"),
        axis.text.x = element_text(face = "bold",size = 14,color="black", vjust=0.5,family = "Calibri"),
        axis.text.y = element_text(face = "bold",size = 14,color="black",family = "Calibri"),
        strip.background = element_blank(),
        legend.position = "none",
        strip.text.x = element_blank(),
        legend.title = element_text(face = "bold",size = 12,color="black",family = "Calibri"),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        legend.text=element_text(face = "bold",size = 12,color="black",family = "Calibri"),
        axis.ticks.length = unit(0.15, "cm"))+
  labs(title = "beta (Lawry [2020] & McCann [2022])",y = "Value", x = "Iteration") +
  annotate(geom = "text",
           x = post_dim(law.mcmc,"saved")/post_dim(law.mcmc,"chains"),
           y = max(MCMCchains(law.mcmc, params = 'beta.law')),
           label = paste("hat(R)","~`=`~",round(as.numeric(fit.law$BUGSoutput$summary[2,8]),3),sep=""),
           hjust = 1.0,vjust=0.0,size = 5.1,family = "Calibri",parse = TRUE)+ #Rhat has to be Changed based on jags output (extraction rounds to nearest interger)
  theme(plot.title = element_text(hjust = 0.5,size = 16,face = "bold",family = "Calibri"))

#### sigma
tracePlot.sigma.law <- ggs_traceplot(law.ggs.dat, family = "sigma.law")+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title.y = element_text(face = "bold", size = 16,vjust = 1,margin = margin(t = 0, r = 50, b = 0, l = 0),family = "Calibri"),
        axis.title.x = element_text(face = "bold", size = 16,vjust = -1,margin = margin(t = 10, r = 0, b = 0, l = 0),family = "Calibri"),
        axis.text.x = element_text(face = "bold",size = 14,color="black", vjust=0.5,family = "Calibri"),
        axis.text.y = element_text(face = "bold",size = 14,color="black",family = "Calibri"),
        strip.background = element_blank(),
        legend.position = "none",
        strip.text.x = element_blank(),
        legend.title = element_text(face = "bold",size = 12,color="black",family = "Calibri"),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        legend.text=element_text(face = "bold",size = 12,color="black",family = "Calibri"),
        axis.ticks.length = unit(0.15, "cm"))+
  labs(title = "sigma (Lawry [2020] & McCann [2022])",y = "Value", x = "Iteration") +
  annotate(geom = "text",
           x = post_dim(law.mcmc,"saved")/post_dim(law.mcmc,"chains"),
           y = max(MCMCchains(law.mcmc, params = 'sigma.law')),
           label = paste("hat(R)","~`=`~",round(as.numeric(fit.law$BUGSoutput$summary[7,8]),3),sep=""),
           hjust = 1.0,vjust=0.0,size = 5.1,family = "Calibri",parse=TRUE)+ # Rhat has to be Changed based on jags output (extraction rounds to nearest interger)
  theme(plot.title = element_text(hjust = 0.5,size = 16,face = "bold",family = "Calibri"))

#### nu
tracePlot.nu.law <- ggs_traceplot(law.ggs.dat, family = "nu.law")+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title.y = element_text(face = "bold", size = 16,vjust = 1,margin = margin(t = 0, r = 50, b = 0, l = 0),family = "Calibri"),
        axis.title.x = element_text(face = "bold", size = 16,vjust = -1,margin = margin(t = 10, r = 0, b = 0, l = 0),family = "Calibri"),
        axis.text.x = element_text(face = "bold",size = 14,color="black", vjust=0.5,family = "Calibri"),
        axis.text.y = element_text(face = "bold",size = 14,color="black",family = "Calibri"),
        strip.background = element_blank(),
        legend.position = "none",
        strip.text.x = element_blank(),
        legend.title = element_text(face = "bold",size = 12,color="black",family = "Calibri"),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        legend.text=element_text(face = "bold",size = 12,color="black",family = "Calibri"),
        axis.ticks.length = unit(0.15, "cm"))+
  labs(title = "nu (Lawry [2020] & McCann [2022])",y = "Value", x = "Iteration") +
  annotate(geom = "text",
           x = post_dim(law.mcmc,"saved")/post_dim(law.mcmc,"chains"),
           y = max(MCMCchains(law.mcmc, params = 'nu.law')),
           label = paste("hat(R)","~`=`~",round(as.numeric(fit.law$BUGSoutput$summary[7,8]),3),sep=""),
           hjust = 1.0,vjust=0.0,size = 5.1,family = "Calibri",parse=TRUE)+  # Rhat has to be Changed based on jags output (extraction rounds to nearest interger)
  theme(plot.title = element_text(hjust = 0.5,size = 16,face = "bold",family = "Calibri"))

### density plots
#### alpha
density.alpha.law <- ggs_density(law.ggs.dat, family = "alpha.law")+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title.y = element_text(face = "bold", size = 16,vjust = 1,margin = margin(t = 0, r = 50, b = 0, l = 0),family = "Calibri"),
        axis.title.x = element_text(face = "bold", size = 16,vjust = -1,margin = margin(t = 10, r = 0, b = 0, l = 0),family = "Calibri"),
        axis.text.x = element_text(face = "bold",size = 14,color="black", vjust=0.5,family = "Calibri"),
        axis.text.y = element_text(face = "bold",size = 14,color="black",family = "Calibri"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.title = element_text(face = "bold",size = 12,color="black",family = "Calibri"),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        legend.text=element_text(face = "bold",size = 12,color="black",family = "Calibri"),
        axis.ticks.length = unit(0.15, "cm"))+
  labs(title = "",y = "Density", x = "Value") +
  theme(plot.title = element_text(hjust = 0.5,size = 16,face = "bold",family = "Calibri"))

#### beta
density.beta.law <- ggs_density(law.ggs.dat, family = "beta.law")+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title.y = element_text(face = "bold", size = 16,vjust = 1,margin = margin(t = 0, r = 50, b = 0, l = 0),family = "Calibri"),
        axis.title.x = element_text(face = "bold", size = 16,vjust = -1,margin = margin(t = 10, r = 0, b = 0, l = 0),family = "Calibri"),
        axis.text.x = element_text(face = "bold",size = 14,color="black", vjust=0.5,family = "Calibri"),
        axis.text.y = element_text(face = "bold",size = 14,color="black",family = "Calibri"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.title = element_text(face = "bold",size = 12,color="black",family = "Calibri"),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        legend.text=element_text(face = "bold",size = 12,color="black",family = "Calibri"),
        axis.ticks.length = unit(0.15, "cm"))+
  labs(title = "",y = "Density", x = "Value") +
  theme(plot.title = element_text(hjust = 0.5,size = 16,face = "bold",family = "Calibri"))

#### sigma
density.sigma.law <- ggs_density(law.ggs.dat, family = "sigma.law")+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title.y = element_text(face = "bold", size = 16,vjust = 1,margin = margin(t = 0, r = 50, b = 0, l = 0),family = "Calibri"),
        axis.title.x = element_text(face = "bold", size = 16,vjust = -1,margin = margin(t = 10, r = 0, b = 0, l = 0),family = "Calibri"),
        axis.text.x = element_text(face = "bold",size = 14,color="black", vjust=0.5,family = "Calibri"),
        axis.text.y = element_text(face = "bold",size = 14,color="black",family = "Calibri"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.title = element_text(face = "bold",size = 12,color="black",family = "Calibri"),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        legend.text=element_text(face = "bold",size = 12,color="black",family = "Calibri"),
        axis.ticks.length = unit(0.15, "cm"))+
  labs(title = "",y = "Density", x = "Value") +
  theme(plot.title = element_text(hjust = 0.5,size = 16,face = "bold",family = "Calibri"))

#### nu
density.nu.law <- ggs_density(law.ggs.dat, family = "nu.law")+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title.y = element_text(face = "bold", size = 16,vjust = 1,margin = margin(t = 0, r = 50, b = 0, l = 0),family = "Calibri"),
        axis.title.x = element_text(face = "bold", size = 16,vjust = -1,margin = margin(t = 10, r = 0, b = 0, l = 0),family = "Calibri"),
        axis.text.x = element_text(face = "bold",size = 14,color="black", vjust=0.5,family = "Calibri"),
        axis.text.y = element_text(face = "bold",size = 14,color="black",family = "Calibri"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.title = element_text(face = "bold",size = 12,color="black",family = "Calibri"),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        legend.text=element_text(face = "bold",size = 12,color="black",family = "Calibri"),
        axis.ticks.length = unit(0.15, "cm"))+
  labs(title = "",y = "Density", x = "Value") +
  theme(plot.title = element_text(hjust = 0.5,size = 16,face = "bold",family = "Calibri"))

### combined plots
#### alpha
combPlot.alphaLaw<-ggarrange(tracePlot.alpha.law,
                          density.alpha.law,
                            ncol = 1,
                            nrow = 2)

png(paste("Output\\Figures\\Bayesian\\Diagnostics\\combPlotLaw.alpha",Sys.Date(),".png",sep=""),
    type="cairo",
    units="in",
    width=9,
    height=12,
    res=300)

print(combPlot.alphaLaw)
dev.off()

#### beta
combPlot.betaLaw<-ggarrange(tracePlot.beta.law,
                          density.beta.law,
                          ncol = 1,
                          nrow = 2)

png(paste("Output\\Figures\\Bayesian\\Diagnostics\\combPlotLaw.beta",Sys.Date(),".png",sep=""),
    type="cairo",
    units="in",
    width=9,
    height=12,
    res=300)

print(combPlot.betaLaw)
dev.off()

#### sigma
combPlot.sigmaLaw<-ggarrange(tracePlot.sigma.law,
                          density.sigma.law,
                          ncol = 1,
                          nrow = 2)

png(paste("Output\\Figures\\Bayesian\\Diagnostics\\combPlotLaw.sigma",Sys.Date(),".png",sep=""),
    type="cairo",
    units="in",
    width=9,
    height=12,
    res=300)

print(combPlot.sigmaLaw)
dev.off()

#### nu
combPlot.nuLaw<-ggarrange(tracePlot.nu.law,
                          density.nu.law,
                          ncol = 1,
                          nrow = 2)

png(paste("Output\\Figures\\Bayesian\\Diagnostics\\combPlotLaw.nu",Sys.Date(),".png",sep=""),
    type="cairo",
    units="in",
    width=9,
    height=12,
    res=300)

print(combPlot.nuLaw)
dev.off()

## generate predictions (fit curve)
mcmc = fit.law$BUGSoutput$sims.matrix

newdata.law = data.frame(xPred.law = seq(0,hiCBP,100))

Xmat.law = model.matrix(~xPred.law, newdata.law)
coefs.law = mcmc[, c("alpha.law", "beta.law")]
fit.predLaw = coefs.law %*% t(Xmat.law)
newdata.law = newdata.law %>% cbind(tidyMCMC(fit.predLaw, conf.int = TRUE, conf.method = "quantile"))

## generate chapter plot
cbpSARequiv.plot_lawry<-ggplot() +
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title.y = element_text(face = "bold", size = 16,vjust = 1,margin = margin(t = 0, r = 10, b = 0, l = 0),family = "Calibri"),
        axis.title.x = element_text(face = "bold", size = 16,vjust = -1,margin = margin(t = 10, r = 0, b = 0, l = 0),family = "Calibri"),
        axis.text.x = element_text(face = "bold",size = 14,color="black", vjust=0.5,family = "Calibri"),
        axis.text.y = element_text(face = "bold",size = 14,color="black",family = "Calibri"),
        legend.title = element_blank(),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        legend.text=element_text(size=12),
        axis.ticks.length = unit(0.15, "cm"))+
  labs(title ="Lawry (2020) & McCann (2022)", y = "SAR", x = "Escapement to Lower Granite Dam") +
  theme(plot.title = element_text(hjust = 0.5,size = 16,face = "bold",family = "Calibri")) +
  geom_ribbon(data =subset(newdata.law,xPred.law >= min(subset(cbpSARequiv.datRed,src=="lawry",select=c(esc))) & xPred.law <= hiCBP),aes(x = xPred.law,ymin = conf.low,ymax = conf.high),alpha=0.3,fill = "#56B4E9")+
  geom_line(data = subset(newdata.law,xPred.law >= min(subset(cbpSARequiv.datRed,src=="lawry",select=c(esc))) & xPred.law <= max(subset(cbpSARequiv.datRed,src=="lawry",select=c(esc)))),aes(xPred.law, estimate), size = 1.0,color = "#56B4E9")+
  geom_line(data = subset(newdata.law,xPred.law > max(subset(cbpSARequiv.datRed,src=="lawry",select=c(esc))) & xPred.law<=hiCBP),aes(xPred.law, estimate), size = 1.0,color = "#56B4E9", linetype = "dashed")+
  geom_segment(aes(x = 0, y = as.numeric(fit.law$BUGSoutput$summary[27,1]), xend = lowCBP), yend = as.numeric(fit.law$BUGSoutput$summary[27,1]), colour = "black",linetype = "dashed",size=1.0)+
  geom_segment(aes(x = lowCBP, y = 0, xend = lowCBP), yend = as.numeric(fit.law$BUGSoutput$summary[27,1]), colour = "black",linetype = "dashed",size=1.0)+
  geom_segment(aes(x = 0, y = as.numeric(fit.law$BUGSoutput$summary[28,1]), xend = medCBP), yend = as.numeric(fit.law$BUGSoutput$summary[28,1]), colour = "black",linetype = "dashed",size=1.0)+
  geom_segment(aes(x = medCBP, y = 0, xend = medCBP), yend = as.numeric(fit.law$BUGSoutput$summary[28,1]), colour = "black",linetype = "dashed",size=1.0)+
  geom_segment(aes(x = 0, y = as.numeric(fit.law$BUGSoutput$summary[26,1]), xend = hiCBP), yend = as.numeric(fit.law$BUGSoutput$summary[26,1]), colour = "black",linetype = "dashed",size=1.0)+
  geom_segment(aes(x = hiCBP, y = 0, xend = hiCBP), yend = as.numeric(fit.law$BUGSoutput$summary[26,1]), colour = "black",linetype = "dashed",size=1.0)+
  geom_point(data = subset(cbpSARequiv.datRed,src=="lawry", select = c(sar,esc)),aes(esc,sar),shape = 21,size = 3.5,stroke=0.5, fill = "#56B4E9")+
  annotate(geom = "text",x = 1000, y = as.numeric(fit.law$BUGSoutput$summary[26,1]), label = "CBP high goal", hjust = 0,vjust=-1,size = 4.5,family = "Calibri")+
  annotate(geom = "text",x = 1000, y = as.numeric(fit.law$BUGSoutput$summary[28,1]), label = "CBP med. goal", hjust = 0,vjust=-1,size = 4.5,family = "Calibri")+
  annotate(geom = "text",x = 1000, y = as.numeric(fit.law$BUGSoutput$summary[27,1]), label = "CBP low goal", hjust = 0,vjust=-1,size = 4.5,family = "Calibri")+
  scale_y_continuous(limits=c(0,0.17),breaks = seq(0,0.17,0.020),labels = function(x) sprintf("%.1f%%", x*100), expand = c(0,0))+
  scale_x_continuous(limits=c(0,hiCBP+10000),breaks = seq(0,hiCBP+10000,25000),labels = function(x) format(x, big.mark = ",", scientific = FALSE), expand = c(0,0))

## output chapter plot
### set image output parameters
png(filename="Output\\Figures\\Bayesian\\Data plots\\cbpSARequivPlot_lawry.png",
    type="cairo",
    units="in",
    width=9,
    height=6,
    res=300)

### write object to working directory
print(cbpSARequiv.plot_lawry)
dev.off()

## define table values
lawLowPred.point <- as.numeric(fit.law$BUGSoutput$summary[27,1])
lawLowPred.lCrI <- as.numeric(fit.law$BUGSoutput$summary[27,3])
lawLowPred.hCrI <- as.numeric(fit.law$BUGSoutput$summary[27,7])
lawMedPred.point <- as.numeric(fit.law$BUGSoutput$summary[28,1])
lawMedPred.lCrI <- as.numeric(fit.law$BUGSoutput$summary[28,3])
lawMedPred.hCrI <- as.numeric(fit.law$BUGSoutput$summary[28,7])
lawHiPred.point <- as.numeric(fit.law$BUGSoutput$summary[26,1])
lawHiPred.lCrI <- as.numeric(fit.law$BUGSoutput$summary[26,3])
lawHiPred.hCrI <- as.numeric(fit.law$BUGSoutput$summary[26,7])

# raymond analysis ----------------------------------------------------------

## define initial and other values
initRay.mod <- lm(data=subset(cbpSARequiv.dat,src=="ray"),sar~esc)
initRay.alpha <- as.numeric(initRay.mod$coef[1])
initRay.alphaSE <- coef(summary(initRay.mod))[1,2]
initRay.beta <- as.numeric(initRay.mod$coef[2])
initRay.betaSE <- coef(summary(initRay.mod))[2,2]
df.ray <- summary(initRay.mod)$df[2]

## Fit model
### formulate data input
cbpSARequiv.rayDat <- list(esc.ray = subset(cbpSARequiv.dat,src=="ray")[,3],
                           sar.ray = subset(cbpSARequiv.dat,src=="ray")[,4],
                           n.ray = nrow(subset(cbpSARequiv.dat,src=="ray")))

### define initial values
inits.ray <- function()
{
  list(alpha.ray = runif(1,initRay.alpha-(5*initRay.alphaSE),initRay.alpha+(5*initRay.alphaSE)),
       beta.ray = runif(1,initRay.beta-(5*initRay.betaSE),initRay.beta+(5*initRay.betaSE)))
}

### specify model
cat('

model {
  # Likelihood
  for (i in 1:n.ray){
    sar.ray[i] ~ dt(mu.ray[i],tau.ray,nu.ray)
    mu.ray[i] <- alpha.ray + beta.ray * esc.ray[i]

    # # res.ray[i] <- sar.ray[i] - mu.ray[i]
    # # sar.new[i] ~ dt(mu.ray[i],tau.ray,nu.ray)
    # # res.sar.new[i] <- sar.new[i] - mu.ray[i]
    #
    # posterior predictive check
    ## predicted values
    ray.pred[i]<-mu.ray[i]
    ## residuals for observed data
    ray.resid[i]<-sar.ray[i]-ray.pred[i]
    ## discrepancy/squared residuals for observed data
    ray.SqResid[i]<-pow(ray.resid[i],2)

    # generate replicate data and compute fit stats
    sar.new.ray[i]~dt(mu.ray[i],tau.ray,nu.ray)
    # discrepancy/squared residuals for new/ideal data
    ray.SqResid.new[i]<-pow(sar.new.ray[i]-ray.pred[i],2)
  }

  # Priors
  alpha.ray ~ dnorm(0 , 0.001)
  beta.ray ~ dnorm(0 , 0.001)
  nuMinusOne.ray ~ dexp(1/29)
  nu.ray <- nuMinusOne.ray+1
  tau.ray <- 1 / (sigma.ray * sigma.ray)
  sigma.ray~dunif(0,100)

  # predict sar at different levels of CBP goals
  ## define CBP levels
  lowCBP <- 43000
  medCBP <- 137000
  hiCBP <- 235000

  ## generate predictions
  muLow.ray <- alpha.ray + beta.ray*lowCBP
  muMed.ray <- alpha.ray + beta.ray*medCBP
  muHi.ray <- alpha.ray + beta.ray*hiCBP
  sarLow.ray ~ dt(muLow.ray,tau.ray,nu.ray)
  sarMed.ray ~ dt(muMed.ray,tau.ray,nu.ray)
  sarHi.ray ~ dt(muHi.ray,tau.ray,nu.ray)

  # derived values
  ## sum of squared residuals for actual dataset
  ray.sumsqO<-sum(ray.SqResid[])
  ## sum of squared residuals for new dataset
  ray.sumsqN<-sum(ray.SqResid.new[])
  ## test whether new data are more extreme
  ray.test<-step(ray.sumsqN-ray.sumsqO)
  ## bayesian pvalue
  ray.bpval<-mean(ray.test)


}', file={lmRay.jags <- tempfile()})

### define parameters to monitor
params.ray <- c("alpha.ray",
                "beta.ray",
                "sigma.ray",
                "nu.ray",
                "sarLow.ray",
                "sarMed.ray",
                "sarHi.ray",
                "ray.sumsqO",
                "ray.sumsqN",
                "ray.bpval",
                "sar.new.ray")

### call jags
fit.ray <- jags(data = cbpSARequiv.rayDat,
                inits = inits.ray,
                parameters.to.save = params.ray,
                model.file = lmRay.jags,
                n.chains = 3,
                n.iter = 250000,
                n.burnin = 25000,
                n.thin = 10,
                DIC = F)

## review output
fit.ray

## review diagnostics
### format output
ray.mcmc <- as.mcmc(fit.ray)
ray.ggs.dat <- ggs(ray.mcmc)

ray.ppc.dat <- data.frame(SumSqO=as.data.frame(subset(ray.ggs.dat,Parameter=="ray.sumsqO",select=c(value))),
                          SumSqN=as.data.frame(subset(ray.ggs.dat,Parameter=="ray.sumsqN",select=c(value))))

colnames(law.ppc.dat) <- c("ray.sumsqO",
                           "ray.sumsqN")
head(ray.ppc.dat)

max(ray.ppc.dat$ray.sumsqN)

### trace plots
#### alpha
tracePlot.alpha.ray <- ggs_traceplot(ray.ggs.dat, family = "alpha.ray")+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title.y = element_text(face = "bold", size = 16,vjust = 1,margin = margin(t = 0, r = 50, b = 0, l = 0),family = "Calibri"),
        axis.title.x = element_text(face = "bold", size = 16,vjust = -1,margin = margin(t = 10, r = 0, b = 0, l = 0),family = "Calibri"),
        axis.text.x = element_text(face = "bold",size = 14,color="black", vjust=0.5,family = "Calibri"),
        axis.text.y = element_text(face = "bold",size = 14,color="black",family = "Calibri"),
        strip.background = element_blank(),
        legend.position = "none",
        strip.text.x = element_blank(),
        legend.title = element_text(face = "bold",size = 12,color="black",family = "Calibri"),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        legend.text=element_text(face = "bold",size = 12,color="black",family = "Calibri"),
        axis.ticks.length = unit(0.15, "cm"))+
  labs(title = "alpha (Raymond [1988])",y = "Value", x = "Iteration") +
  annotate(geom = "text",
           x = post_dim(ray.mcmc,"saved")/post_dim(ray.mcmc,"chains"),
           y = max(MCMCchains(ray.mcmc, params = 'alpha.ray')),
           label = paste("hat(R)","~`=`~",round(as.numeric(fit.ray$BUGSoutput$summary[1,8]),3),sep=""),
           hjust = 1.0,vjust=0.0,size = 5.1,family = "Calibri",parse = TRUE)+
  theme(plot.title = element_text(hjust = 0.5,size = 16,face = "bold",family = "Calibri"))

#### beta
tracePlot.beta.ray <- ggs_traceplot(ray.ggs.dat, family = "beta.ray")+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title.y = element_text(face = "bold", size = 16,vjust = 1,margin = margin(t = 0, r = 50, b = 0, l = 0),family = "Calibri"),
        axis.title.x = element_text(face = "bold", size = 16,vjust = -1,margin = margin(t = 10, r = 0, b = 0, l = 0),family = "Calibri"),
        axis.text.x = element_text(face = "bold",size = 14,color="black", vjust=0.5,family = "Calibri"),
        axis.text.y = element_text(face = "bold",size = 14,color="black",family = "Calibri"),
        strip.background = element_blank(),
        legend.position = "none",
        strip.text.x = element_blank(),
        legend.title = element_text(face = "bold",size = 12,color="black",family = "Calibri"),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        legend.text=element_text(face = "bold",size = 12,color="black",family = "Calibri"),
        axis.ticks.length = unit(0.15, "cm"))+
  labs(title = "beta (Raymond [1988])",y = "Value", x = "Iteration") +
  annotate(geom = "text",
           x = post_dim(ray.mcmc,"saved")/post_dim(ray.mcmc,"chains"),
           y = max(MCMCchains(ray.mcmc, params = 'beta.ray')),
           label = paste("hat(R)","~`=`~",round(as.numeric(fit.ray$BUGSoutput$summary[2,8]),3),sep=""),
           hjust = 1.0,vjust=0.0,size = 5.1,family = "Calibri",parse = TRUE)+
  theme(plot.title = element_text(hjust = 0.5,size = 16,face = "bold",family = "Calibri"))

#### sigma
tracePlot.sigma.ray <- ggs_traceplot(ray.ggs.dat, family = "sigma.ray")+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title.y = element_text(face = "bold", size = 16,vjust = 1,margin = margin(t = 0, r = 50, b = 0, l = 0),family = "Calibri"),
        axis.title.x = element_text(face = "bold", size = 16,vjust = -1,margin = margin(t = 10, r = 0, b = 0, l = 0),family = "Calibri"),
        axis.text.x = element_text(face = "bold",size = 14,color="black", vjust=0.5,family = "Calibri"),
        axis.text.y = element_text(face = "bold",size = 14,color="black",family = "Calibri"),
        strip.background = element_blank(),
        legend.position = "none",
        strip.text.x = element_blank(),
        legend.title = element_text(face = "bold",size = 12,color="black",family = "Calibri"),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        legend.text=element_text(face = "bold",size = 12,color="black",family = "Calibri"),
        axis.ticks.length = unit(0.15, "cm"))+
  labs(title = "sigma (Raymond [1988])",y = "Value", x = "Iteration") +
  annotate(geom = "text",
           x = post_dim(ray.mcmc,"saved")/post_dim(ray.mcmc,"chains"),
           y = max(MCMCchains(ray.mcmc, params = 'sigma.ray')),
           label = paste("hat(R)","~`=`~",round(as.numeric(fit.ray$BUGSoutput$summary[7,8]),3),sep=""),
           hjust = 1.0,vjust=0.0,size = 5.1,family = "Calibri",parse=TRUE)+
  theme(plot.title = element_text(hjust = 0.5,size = 16,face = "bold",family = "Calibri"))

#### nu
tracePlot.nu.ray <- ggs_traceplot(ray.ggs.dat, family = "nu.ray")+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title.y = element_text(face = "bold", size = 16,vjust = 1,margin = margin(t = 0, r = 50, b = 0, l = 0),family = "Calibri"),
        axis.title.x = element_text(face = "bold", size = 16,vjust = -1,margin = margin(t = 10, r = 0, b = 0, l = 0),family = "Calibri"),
        axis.text.x = element_text(face = "bold",size = 14,color="black", vjust=0.5,family = "Calibri"),
        axis.text.y = element_text(face = "bold",size = 14,color="black",family = "Calibri"),
        strip.background = element_blank(),
        legend.position = "none",
        strip.text.x = element_blank(),
        legend.title = element_text(face = "bold",size = 12,color="black",family = "Calibri"),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        legend.text=element_text(face = "bold",size = 12,color="black",family = "Calibri"),
        axis.ticks.length = unit(0.15, "cm"))+
  labs(title = "nu (Raymond [1988])",y = "Value", x = "Iteration") +
  annotate(geom = "text",
           x = post_dim(ray.mcmc,"saved")/post_dim(ray.mcmc,"chains"),
           y = max(MCMCchains(ray.mcmc, params = 'nu.ray')),
           label = paste("hat(R)","~`=`~",round(as.numeric(fit.ray$BUGSoutput$summary[7,8]),3),sep=""),
           hjust = 1.0,vjust=0.0,size = 5.1,family = "Calibri",parse=TRUE)+
  theme(plot.title = element_text(hjust = 0.5,size = 16,face = "bold",family = "Calibri"))

### density plots
#### alpha
density.alpha.ray <- ggs_density(ray.ggs.dat, family = "alpha.ray")+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title.y = element_text(face = "bold", size = 16,vjust = 1,margin = margin(t = 0, r = 50, b = 0, l = 0),family = "Calibri"),
        axis.title.x = element_text(face = "bold", size = 16,vjust = -1,margin = margin(t = 10, r = 0, b = 0, l = 0),family = "Calibri"),
        axis.text.x = element_text(face = "bold",size = 14,color="black", vjust=0.5,family = "Calibri"),
        axis.text.y = element_text(face = "bold",size = 14,color="black",family = "Calibri"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.title = element_text(face = "bold",size = 12,color="black",family = "Calibri"),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        legend.text=element_text(face = "bold",size = 12,color="black",family = "Calibri"),
        axis.ticks.length = unit(0.15, "cm"))+
  labs(title = "",y = "Density", x = "Value") +
  theme(plot.title = element_text(hjust = 0.5,size = 16,face = "bold",family = "Calibri"))

#### beta
density.beta.ray <- ggs_density(ray.ggs.dat, family = "beta.ray")+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title.y = element_text(face = "bold", size = 16,vjust = 1,margin = margin(t = 0, r = 50, b = 0, l = 0),family = "Calibri"),
        axis.title.x = element_text(face = "bold", size = 16,vjust = -1,margin = margin(t = 10, r = 0, b = 0, l = 0),family = "Calibri"),
        axis.text.x = element_text(face = "bold",size = 14,color="black", vjust=0.5,family = "Calibri"),
        axis.text.y = element_text(face = "bold",size = 14,color="black",family = "Calibri"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.title = element_text(face = "bold",size = 12,color="black",family = "Calibri"),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        legend.text=element_text(face = "bold",size = 12,color="black",family = "Calibri"),
        axis.ticks.length = unit(0.15, "cm"))+
  labs(title = "",y = "Density", x = "Value") +
  theme(plot.title = element_text(hjust = 0.5,size = 16,face = "bold",family = "Calibri"))

#### sigma
density.sigma.ray <- ggs_density(ray.ggs.dat, family = "sigma.ray")+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title.y = element_text(face = "bold", size = 16,vjust = 1,margin = margin(t = 0, r = 50, b = 0, l = 0),family = "Calibri"),
        axis.title.x = element_text(face = "bold", size = 16,vjust = -1,margin = margin(t = 10, r = 0, b = 0, l = 0),family = "Calibri"),
        axis.text.x = element_text(face = "bold",size = 14,color="black", vjust=0.5,family = "Calibri"),
        axis.text.y = element_text(face = "bold",size = 14,color="black",family = "Calibri"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.title = element_text(face = "bold",size = 12,color="black",family = "Calibri"),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        legend.text=element_text(face = "bold",size = 12,color="black",family = "Calibri"),
        axis.ticks.length = unit(0.15, "cm"))+
  labs(title = "",y = "Density", x = "Value") +
  theme(plot.title = element_text(hjust = 0.5,size = 16,face = "bold",family = "Calibri"))

#### nu
density.nu.ray <- ggs_density(ray.ggs.dat, family = "nu.ray")+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title.y = element_text(face = "bold", size = 16,vjust = 1,margin = margin(t = 0, r = 50, b = 0, l = 0),family = "Calibri"),
        axis.title.x = element_text(face = "bold", size = 16,vjust = -1,margin = margin(t = 10, r = 0, b = 0, l = 0),family = "Calibri"),
        axis.text.x = element_text(face = "bold",size = 14,color="black", vjust=0.5,family = "Calibri"),
        axis.text.y = element_text(face = "bold",size = 14,color="black",family = "Calibri"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.title = element_text(face = "bold",size = 12,color="black",family = "Calibri"),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        legend.text=element_text(face = "bold",size = 12,color="black",family = "Calibri"),
        axis.ticks.length = unit(0.15, "cm"))+
  labs(title = "",y = "Density", x = "Value") +
  theme(plot.title = element_text(hjust = 0.5,size = 16,face = "bold",family = "Calibri"))

### combined plots
#### alpha
combPlot.alphaRay<-ggarrange(tracePlot.alpha.ray,
                             density.alpha.ray,
                             ncol = 1,
                             nrow = 2)

png(paste("Output\\Figures\\Bayesian\\Diagnostics\\combPlotRay.alpha",Sys.Date(),".png",sep=""),
    type="cairo",
    units="in",
    width=9,
    height=12,
    res=300)

print(combPlot.alphaRay)
dev.off()

#### beta
combPlot.betaRay<-ggarrange(tracePlot.beta.ray,
                            density.beta.ray,
                            ncol = 1,
                            nrow = 2)

png(paste("Output\\Figures\\Bayesian\\Diagnostics\\combPlotRay.beta",Sys.Date(),".png",sep=""),
    type="cairo",
    units="in",
    width=9,
    height=12,
    res=300)

print(combPlot.betaRay)
dev.off()

#### sigma
combPlot.sigmaRay<-ggarrange(tracePlot.sigma.ray,
                             density.sigma.ray,
                             ncol = 1,
                             nrow = 2)

png(paste("Output\\Figures\\Bayesian\\Diagnostics\\combPlotRay.sigma",Sys.Date(),".png",sep=""),
    type="cairo",
    units="in",
    width=9,
    height=12,
    res=300)

print(combPlot.sigmaRay)
dev.off()

#### nu
combPlot.nuRay<-ggarrange(tracePlot.nu.ray,
                          density.nu.ray,
                          ncol = 1,
                          nrow = 2)

png(paste("Output\\Figures\\Bayesian\\Diagnostics\\combPlotRay.nu",Sys.Date(),".png",sep=""),
    type="cairo",
    units="in",
    width=9,
    height=12,
    res=300)

print(combPlot.nuRay)
dev.off()

## generate predictions (fit curve)
mcmc.ray = fit.ray$BUGSoutput$sims.matrix

newdata.ray = data.frame(xPred.ray = seq(0,hiCBP,100))

Xmat.ray = model.matrix(~xPred.ray, newdata.ray)
coefs.ray = mcmc.ray[, c("alpha.ray", "beta.ray")]
fit.predRay = coefs.ray %*% t(Xmat.ray)
newdata.ray = newdata.ray %>% cbind(tidyMCMC(fit.predRay, conf.int = TRUE, conf.method = "quantile"))

## generate chapter plot
cbpSARequiv.plot_raymond<-ggplot() +
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title.y = element_text(face = "bold", size = 16,vjust = 1,margin = margin(t = 0, r = 10, b = 0, l = 0),family = "Calibri"),
        axis.title.x = element_text(face = "bold", size = 16,vjust = -1,margin = margin(t = 10, r = 0, b = 0, l = 0),family = "Calibri"),
        axis.text.x = element_text(face = "bold",size = 14,color="black", vjust=0.5,family = "Calibri"),
        axis.text.y = element_text(face = "bold",size = 14,color="black",family = "Calibri"),
        legend.title = element_blank(),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        legend.text=element_text(size=12),
        axis.ticks.length = unit(0.15, "cm"))+
  labs(title ="Raymond (1998)", y = "SAR", x = "") +
  theme(plot.title = element_text(hjust = 0.5,size = 16,face = "bold",family = "Calibri")) +
  geom_ribbon(data =subset(newdata.ray,xPred.ray >= min(subset(cbpSARequiv.datRed,src=="ray",select=c(esc))) & xPred.ray <= hiCBP),aes(x = xPred.ray,ymin = conf.low,ymax = conf.high),alpha=0.3,fill = "#E69F00")+
  geom_line(data = subset(newdata.ray,xPred.ray >= min(subset(cbpSARequiv.datRed,src=="ray",select=c(esc))) & xPred.ray <= max(subset(cbpSARequiv.datRed,src=="ray",select=c(esc)))),aes(xPred.ray, estimate), size = 1.0,color = "#E69F00")+
  geom_line(data = subset(newdata.ray,xPred.ray > max(subset(cbpSARequiv.datRed,src=="ray",select=c(esc))) & xPred.ray<=hiCBP),aes(xPred.ray, estimate), size = 1.0,color = "#E69F00", linetype = "dashed")+
  geom_segment(aes(x = 0, y = as.numeric(fit.ray$BUGSoutput$summary[29,1]), xend = lowCBP), yend = as.numeric(fit.ray$BUGSoutput$summary[29,1]), colour = "black",linetype = "dashed",size=1.0)+
  geom_segment(aes(x = lowCBP, y = 0, xend = lowCBP), yend = as.numeric(fit.ray$BUGSoutput$summary[29,1]), colour = "black",linetype = "dashed",size=1.0)+
  geom_segment(aes(x = 0, y = as.numeric(fit.ray$BUGSoutput$summary[30,1]), xend = medCBP), yend = as.numeric(fit.ray$BUGSoutput$summary[30,1]), colour = "black",linetype = "dashed",size=1.0)+
  geom_segment(aes(x = medCBP, y = 0, xend = medCBP), yend = as.numeric(fit.ray$BUGSoutput$summary[30,1]), colour = "black",linetype = "dashed",size=1.0)+
  geom_segment(aes(x = 0, y = as.numeric(fit.ray$BUGSoutput$summary[28,1]), xend = hiCBP), yend = as.numeric(fit.ray$BUGSoutput$summary[28,1]), colour = "black",linetype = "dashed",size=1.0)+
  geom_segment(aes(x = hiCBP, y = 0, xend = hiCBP), yend = as.numeric(fit.ray$BUGSoutput$summary[28,1]), colour = "black",linetype = "dashed",size=1.0)+
  geom_point(data = subset(cbpSARequiv.datRed,src=="ray", select = c(sar,esc)),aes(esc,sar),shape = 24,size = 3.5,stroke=0.5, fill = "#E69F00")+
  annotate(geom = "text",x = 1000, y = as.numeric(fit.ray$BUGSoutput$summary[28,1]), label = "CBP high goal", hjust = 0,vjust=-1,size = 4.5,family = "Calibri")+
  annotate(geom = "text",x = 1000, y = as.numeric(fit.ray$BUGSoutput$summary[30,1]), label = "CBP med. goal", hjust = 0,vjust=-1,size = 4.5,family = "Calibri")+
  annotate(geom = "text",x = 1000, y = as.numeric(fit.ray$BUGSoutput$summary[29,1]), label = "CBP low goal", hjust = 0,vjust=-1,size = 4.5,family = "Calibri")+
  scale_y_continuous(limits=c(0,0.17),breaks = seq(0,0.17,0.020),labels = function(x) sprintf("%.1f%%", x*100), expand = c(0,0))+
  scale_x_continuous(limits=c(0,hiCBP+10000),breaks = seq(0,hiCBP+10000,25000),labels = function(x) format(x, big.mark = ",", scientific = FALSE), expand = c(0,0))

## output chapter plot
### set image output parameters
png(filename="Output\\Figures\\Bayesian\\Data plots\\cbpSARequivPlot_raymond.png",
    type="cairo",
    units="in",
    width=9,
    height=6,
    res=300)

### write object to working directory
print(cbpSARequiv.plot_raymond)
dev.off()

## define table values
rayLowPred.point <- as.numeric(fit.ray$BUGSoutput$summary[29,1])
rayLowPred.lCrI <- as.numeric(fit.ray$BUGSoutput$summary[29,3])
rayLowPred.hCrI <- as.numeric(fit.ray$BUGSoutput$summary[29,7])
rayMedPred.point <- as.numeric(fit.ray$BUGSoutput$summary[30,1])
rayMedPred.lCrI <- as.numeric(fit.ray$BUGSoutput$summary[30,3])
rayMedPred.hCrI <- as.numeric(fit.ray$BUGSoutput$summary[30,7])
rayHiPred.point <- as.numeric(fit.ray$BUGSoutput$summary[28,1])
rayHiPred.lCrI <- as.numeric(fit.ray$BUGSoutput$summary[28,3])
rayHiPred.hCrI <- as.numeric(fit.ray$BUGSoutput$summary[28,7])

#### Combined figure
cbpSARequiv.plot_comb<-ggarrange(cbpSARequiv.plot_raymond,
                                 cbpSARequiv.plot_lawry,
                                 ncol = 1,
                                 nrow = 2)

png(paste("Output\\Figures\\Bayesian\\Data plots\\cbpSARequivPlot_comb",Sys.Date(),".png",sep=""),
    type="cairo",
    units="in",
    width=9,
    height=12,
    res=300)

print(cbpSARequiv.plot_comb)
dev.off()