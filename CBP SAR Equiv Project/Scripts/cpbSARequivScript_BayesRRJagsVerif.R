# call packages -----------------------------------------------------------
packages <- c("snow",
              "snowfall",
              "coda",
              "HDInterval",
              "mcmcplots",
              "car",
              "openair",
              "ggplot2",
              "zoo",
              "openxlsx")

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
cbpSARequiv.verifDat <- read.xlsx("Input\\SARvLWGsp-su_esc.xlsx",
                             sheet = 1,
                             colNames = TRUE)

# lawry analysis ----------------------------------------------------------
lawPredPerfFun.RunOut<-data.frame(Pred=numeric(),
                                  Obs=numeric(),
                                  stringsAsFactors = FALSE)

lawPredPerfFun <- function(nRec){
  ### formulate data input
  cbpSARequiv.lawDatVerif <- list(esc.law.verif = c(subset(cbpSARequiv.verifDat,src=="lawry")[-nRec,3],subset(cbpSARequiv.verifDat,src=="lawry")[nRec,3]),
                             sar.law.verif = c(subset(cbpSARequiv.verifDat,src=="lawry")[-nRec,4],NA),
                             n.law.verif = nrow(subset(cbpSARequiv.verifDat,src=="lawry")))

## define initial and other values
initLaw.modVerif <- lm(data=subset(cbpSARequiv.verifDat,src=="lawry")[-nRec,],sar~esc)
initLaw.alphaVerif <- as.numeric(initLaw.modVerif$coef[1])
initLaw.alphaSEverif <- coef(summary(initLaw.modVerif))[1,2]
initLaw.betaVerif <- as.numeric(initLaw.modVerif$coef[2])
initLaw.betaSEverif <- coef(summary(initLaw.modVerif))[2,2]
df.lawVerrif <- summary(initLaw.modVerif)$df[2]

### define initial values
inits.law <- function()
{
  list(alpha.law = runif(1,initLaw.alphaVerif-(5*initLaw.alphaSEverif),initLaw.alphaVerif+(5*initLaw.alphaSEverif)),
       beta.law = runif(1,initLaw.betaVerif-(5*initLaw.betaSEverif),initLaw.betaVerif+(5*initLaw.betaSEverif)))
}

### specify model
cat('

model {
  # Likelihood
  for (i in 1:n.law.verif){
    sar.law.verif[i] ~ dt(mu.law.verif[i],tau.law.verif,nu.law.verif)
    mu.law.verif[i] <- alpha.law.verif + beta.law.verif * esc.law.verif[i]
  }

  # Priors
  alpha.law.verif ~ dnorm(0 , 0.001)
  beta.law.verif ~ dnorm(0 , 0.001)
  nuMinusOne.law.verif ~ dexp(1/29)
  nu.law.verif <- nuMinusOne.law.verif+1
  tau.law.verif <- 1 / (sigma.law.verif * sigma.law.verif)
  sigma.law.verif~dunif(0,100)

}', file={lmLaw.jags.verif <- tempfile()})

### define parameters to monitor
params.law <- c("mu.law.verif")

### call jags
fit.law.verif <- jags(data = cbpSARequiv.lawDatVerif,
                inits = inits.law,
                parameters.to.save = params.law,
                model.file = lmLaw.jags.verif,
                n.chains = 3,
                n.iter = 250000,
                n.burnin = 25000,
                n.thin = 10,
                DIC = F)

# ## review output
# fit.law

fit.law.mcmc.verif <- as.mcmc(fit.law.verif)
# combine output ----------------------------------------------------------

lawPredPerfFun.RunOut[nrow(lawPredPerfFun.RunOut)+1,]<<-c(as.numeric(summary(fit.law.mcmc.verif)$statistics[nrow(subset(cbpSARequiv.verifDat,src=="lawry")),1]),
                                                          subset(cbpSARequiv.verifDat,src=="lawry")[nRec,4])

}

for(nRec in seq(1,nrow(subset(cbpSARequiv.verifDat,src=="lawry")),1)){
  lawPredPerfFun(nRec)
}

IOA.law2 <- modStats(lawPredPerfFun.RunOut,mod = "Pred",obs = "Obs",statistic = c("IOA"))

cbpSARequiv.PPplot_lawry<-ggplot() +
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title.y = element_text(face = "bold", size = 16,vjust = 1,margin = margin(t = 0, r = 10, b = 0, l = 0),family = "Calibri"),
        axis.title.x = element_text(face = "bold", size = 16,vjust = -1,margin = margin(t = 10, r = 0, b = 0, l = 0),family = "Calibri"),
        axis.text.x = element_text(face = "bold",size = 16,color="black", vjust=0.5,family = "Calibri"),
        axis.text.y = element_text(face = "bold",size = 16,color="black",family = "Calibri"),
        legend.title = element_blank(),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        legend.text=element_text(size=12),
        axis.ticks.length = unit(0.15, "cm"))+
  labs(title ="Lawry (2020) & McCann (2022)", y = "Estimated", x = "Mean posterior predictions") +
  theme(plot.title = element_text(hjust = 0.5,size = 16,face = "bold",family = "Calibri")) +
  geom_abline(intercept = 0, slope = 1)+
  geom_point(data = lawPredPerfFun.RunOut,aes(Pred,Obs),shape = 21,size = 3.5,stroke=0.5, fill = "#56B4E9")+
  scale_y_continuous(limits=c(0,0.05),breaks = seq(0,0.05,0.01),labels = scales::percent, expand = c(0,0))+
  scale_x_continuous(limits=c(0,0.05),breaks = seq(0,0.05,0.01),labels = scales::percent, expand = c(0,0))+
  annotate(geom = "text",x = 0.001, y = 0.047, label = bquote(italic(d[r]~"="~.(round(as.numeric(IOA.law2[,2]),2)))),hjust=0,size = 5.1,family = "Calibri")


png(filename=paste("Output\\Figures\\Bayesian\\Predictive performance\\cbpSARequivPPplot_lawryBayes",Sys.Date(),".png",sep=""),
    type="cairo",
    units="in",
    width=8,
    height=6,
    res=300)

print(cbpSARequiv.PPplot_lawry)
dev.off() #  turn device off

# raymond analysis ----------------------------------------------------------

rayPredPerfFun.RunOut<-data.frame(Pred=numeric(),
                                  Obs=numeric(),
                                  stringsAsFactors = FALSE)

rayPredPerfFun <- function(nRec){
  ### formulate data input
  cbpSARequiv.rayDatVerif <- list(esc.ray.verif = c(subset(cbpSARequiv.verifDat,src=="ray")[-nRec,3],subset(cbpSARequiv.verifDat,src=="ray")[nRec,3]),
                                  sar.ray.verif = c(subset(cbpSARequiv.verifDat,src=="ray")[-nRec,4],NA),
                                  n.ray.verif = nrow(subset(cbpSARequiv.verifDat,src=="ray")))

  ## define initial and other values
  initray.modVerif <- lm(data=subset(cbpSARequiv.verifDat,src=="ray")[-nRec,],sar~esc)
  initray.alphaVerif <- as.numeric(initray.modVerif$coef[1])
  initray.alphaSEverif <- coef(summary(initray.modVerif))[1,2]
  initray.betaVerif <- as.numeric(initray.modVerif$coef[2])
  initray.betaSEverif <- coef(summary(initray.modVerif))[2,2]
  df.rayVerrif <- summary(initray.modVerif)$df[2]

  ### define initial values
  inits.ray <- function()
  {
    list(alpha.ray = runif(1,initray.alphaVerif-(5*initray.alphaSEverif),initray.alphaVerif+(5*initray.alphaSEverif)),
         beta.ray = runif(1,initray.betaVerif-(5*initray.betaSEverif),initray.betaVerif+(5*initray.betaSEverif)))
  }

  ### specify model
  cat('

model {
  # Likelihood
  for (i in 1:n.ray.verif){
    sar.ray.verif[i] ~ dt(mu.ray.verif[i],tau.ray.verif,nu.ray.verif)
    mu.ray.verif[i] <- alpha.ray.verif + beta.ray.verif * esc.ray.verif[i]
  }

  # Priors
  alpha.ray.verif ~ dnorm(0 , 0.001)
  beta.ray.verif ~ dnorm(0 , 0.001)
  nuMinusOne.ray.verif ~ dexp(1/29)
  nu.ray.verif <- nuMinusOne.ray.verif+1
  tau.ray.verif <- 1 / (sigma.ray.verif * sigma.ray.verif)
  sigma.ray.verif~dunif(0,100)

}', file={lmray.jags.verif <- tempfile()})

  ### define parameters to monitor
  params.ray <- c("mu.ray.verif")

  ### call jags
  fit.ray.verif <- jags(data = cbpSARequiv.rayDatVerif,
                        inits = inits.ray,
                        parameters.to.save = params.ray,
                        model.file = lmray.jags.verif,
                        n.chains = 3,
                        n.iter = 250000,
                        n.burnin = 25000,
                        n.thin = 10,
                        DIC = F)

  # ## review output
  # fit.ray

  fit.ray.mcmc.verif <- as.mcmc(fit.ray.verif)
  # combine output ----------------------------------------------------------

  rayPredPerfFun.RunOut[nrow(rayPredPerfFun.RunOut)+1,]<<-c(as.numeric(summary(fit.ray.mcmc.verif)$statistics[nrow(subset(cbpSARequiv.verifDat,src=="ray")),1]),
                                                            subset(cbpSARequiv.verifDat,src=="ray")[nRec,4])

}

for(nRec in seq(1,nrow(subset(cbpSARequiv.verifDat,src=="ray")),1)){
  rayPredPerfFun(nRec)
}

IOA.ray2 <- modStats(rayPredPerfFun.RunOut,mod = "Pred",obs = "Obs",statistic = c("IOA"))

cbpSARequiv.PPplot_raymond<-ggplot() +
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title.y = element_text(face = "bold", size = 16,vjust = 1,margin = margin(t = 0, r = 10, b = 0, l = 0),family = "Calibri"),
        axis.title.x = element_text(face = "bold", size = 16,vjust = -1,margin = margin(t = 10, r = 0, b = 0, l = 0),family = "Calibri"),
        axis.text.x = element_text(face = "bold",size = 16,color="black", vjust=0.5,family = "Calibri"),
        axis.text.y = element_text(face = "bold",size = 16,color="black",family = "Calibri"),
        legend.title = element_blank(),
        plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"),
        legend.text=element_text(size=12),
        axis.ticks.length = unit(0.15, "cm"))+
  labs(title ="Raymond (1988)", y = "Estimated", x = "Mean posterior predictions") +
  theme(plot.title = element_text(hjust = 0.5,size = 16,face = "bold",family = "Calibri")) +
  geom_abline(intercept = 0, slope = 1)+
  geom_point(data = rayPredPerfFun.RunOut,aes(Pred,Obs),shape = 24,size = 3.5,stroke=0.5, fill = "#E69F00")+
  scale_y_continuous(limits=c(0,0.05),breaks = seq(0,0.06,0.01),labels = scales::percent, expand = c(0,0))+
  scale_x_continuous(limits=c(0,0.05),breaks = seq(0,0.05,0.01),labels = scales::percent, expand = c(0,0))+
  annotate(geom = "text",x = 0.001, y = 0.047, label = bquote(italic(d[r]~"="~.(round(as.numeric(IOA.ray2[,2]),2)))),hjust=0,size = 5.1,family = "Calibri")

png(filename=paste("Output\\Figures\\Bayesian\\Predictive performance\\cbpSARequivPPplot_raymondBayes",Sys.Date(),".png",sep=""),
    type="cairo",
    units="in",
    width=8,
    height=6,
    res=300)

print(cbpSARequiv.PPplot_raymond)
dev.off() #  turn device off

# combined figures --------------------------------------------------------
cbpSARequiv.plotPP_comb<-ggarrange(cbpSARequiv.PPplot_raymond,
                                   cbpSARequiv.PPplot_lawry,
                                   ncol = 1,
                                   nrow = 2)

png(paste("Output\\Figures\\Bayesian\\Predictive performance\\cbpSARequivPlotPP_combBayes",Sys.Date(),".png",sep=""),
    type="cairo",
    units="in",
    width=9,
    height=12,
    res=300)

print(cbpSARequiv.plotPP_comb)
dev.off()



