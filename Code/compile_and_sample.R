#This file runs the code in the file "land_change_analysis_v2.R" to load the data
#and three models to be fit to it, then compiles and samples the models.

#load the data and models
source("Code/land_change_analysis_v2.R")
#load packages
library(rjags)

####compile models and take samples using 10 chains ####
##model 1
mod_1<-
  jags.model(
    textConnection(no_pool),
    data = dat_mod1_mod2,
    n.chains = 10,
    n.adapt = 1000
  )
update(mod_1, 5000)
mod1_dic <- dic.samples(model = mod_1, n.iter = 5000, type = "popt")

##model 2
mod_2<-
  jags.model(
    textConnection(coast_pool),
    data = dat_mod1_mod2,
    n.chains = 10,
    n.adapt = 1000
  )
update(mod_2, 5000)
#sample dic
mod2_dic <- dic.samples(model = mod_2, n.iter = 5000, type = "popt")

##model 3
mod_3<-
  jags.model(
    textConnection(province_coast_pool),
    data = dat_mod3,
    n.chains = 10,
    n.adapt = 1000
  )
update(mod_3, 5000)
#sample the site-level 
  mod3_samples_site<-coda.samples(mod_3, c("bs", "bi", "bh","bx"), n.iter = 30000,thin = 5)
  saveRDS(mod3_samples_site,file="Output/Results/mod3_site_samples_beta.rds")
#sample the site-level intercepts (alpha)
  mod3_samples_site_alpha<-coda.samples(mod_3, c("alpha"), n.iter = 30000,thin = 5)
  saveRDS(mod3_samples_site_alpha,file="Output/Results/mod3_site_samples_intercepts.rds")
#sample at the prov scale
  mod3_samples_prov<-coda.samples(mod_3, c("b.prov.h","b.prov.i","b.prov.s","b.prov.x"), n.iter = 30000,thin = 5)
  saveRDS(mod3_samples_prov,file="Output/Results/mod3_prov_samples.rds")

#sample predicted values
  mod3_samples_mu<-coda.samples(mod_3, "mu", n.iter = 30000,thin = 5)
  saveRDS(mod3_samples_mu,file="Output/Results/mod3_mu_samples.rds")
#get points for the interaction plot
  #mod3_samples_interaction_plot<-coda.samples(mod_3, c("c.low", "c.high", "d.low","d.high"), n.iter = 10000,thin = 5)
  mod3_samples_interaction_plot<-coda.samples(mod_3, c("c.low2", "c.high2", "d.low2","d.high2"), n.iter = 10000,thin = 5)
  saveRDS(mod3_samples_interaction_plot,file="Output/Results/mod3_interaction_plot.rds") 
#sample dic
mod3_dic <- dic.samples(model = mod_3, n.iter = 5000, type = "popt")


###compare models with dic (lower is better) ###
mod1_dic
mod2_dic
mod3_dic

diffdic(mod1_dic,mod2_dic)
diffdic(mod1_dic,mod3_dic)
diffdic(mod2_dic,mod3_dic)


###look at model 2 site level estimates ###
z<-summary(mod2_samples_site)
which(z$quantiles[,"2.5%"]<0&z$quantiles[,"97.5%"]<0)
which(z$quantiles[,"2.5%"]>0&z$quantiles[,"97.5%"]>0)

z3<-summary(mod3_samples_site)
which(z3$quantiles[,"2.5%"]<0&z3$quantiles[,"97.5%"]<0)
which(z3$quantiles[,"2.5%"]>0&z3$quantiles[,"97.5%"]>0)

x11()
plot(mod2_samples_site,ask = T)
z_mu<-summary(mod2_samples_mu)$quantiles[,"50%"]
plot(df1$land_change,z_mu)
cor(df1$land_change,z_mu)

#test convergence of site-level parameters
len<-dim(mod3_samples_site_alpha[[1]])[2]
N<-dim(mod3_samples_site_alpha[[1]])[1]
M<-length(mod3_samples_site_alpha)
G_R_conv<-rep(NA,len)
for(i in 1:len){
  tmp<-sapply(1:M,function(x)mod3_samples_site_alpha[[x]][,i])
  B<-N/(M-1)*var(apply(tmp,2,mean))
  apply(tmp,2,var)
  W<-mean(apply(tmp,2,var))
  V<-(N-1)/N*W+(M+1)/(M*N)*B
  G_R_conv[i]<-sqrt(V/W)
}
#these should all be around 1
quantile(G_R_conv)

#coast-level convergence
len<-dim(mod2_samples_coast[[1]])[2]
N<-dim(mod2_samples_coast[[1]])[1]
M<-length(mod2_samples_coast)
G_R_conv<-rep(NA,len)
for(i in 1:len){
  tmp<-sapply(1:M,function(x)mod2_samples_coast[[x]][,i])
  B<-N/(M-1)*var(apply(tmp,2,mean))
  apply(tmp,2,var)
  W<-mean(apply(tmp,2,var))
  V<-(N-1)/N*W+(M+1)/(M*N)*B
  G_R_conv[i]<-sqrt(V/W)
}
quantile(G_R_conv)


