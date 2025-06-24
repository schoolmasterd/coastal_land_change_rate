###This file contains the models and data for the analysis###
### add title etc ###

#load data
df1 <- read.csv("Data/combined_data_for_analysis.csv")


#set up one-hot coding for basin-level analysis
bas <- matrix(0, nrow = dim(df1)[1], ncol = 9)
colnames(bas) <- unique(df1$basin)
for (i in 1:dim(df1)[1])
  bas[i, df1$basin[i]] <- 1
#make site-level lookup table for basin
blu <- 1:9
names(blu) <- unique(df1$basin)
bas_by_site <- blu[df1$basin[!duplicated(df1$site)]]
prov_by_site<- prov_lookup[df1$basin[!duplicated(df1$site)]]
saveRDS(prov_by_site,"Output/Results/prov_by_site.RData")


saveRDS(bas_by_site,file = "Output/Results/basin_by_site.rds")
site<-unique(df1$site)
saveRDS(site,file = "Output/Results/site_list.rds")


#Here we set up three (JAGS) models of different hierarchical structure for comparison
#one will have not pooling, i.e., each site will be fit independently of all the others
#in the second, we will have all the site-level estimates drawn from coastwide-level
#distribution. Finally, we will draw site-level estimates from basin-scale distributions
#and the basin scale parameters from a coastwide-scale distribution.
#we will use DIC to compare the fits of the models
#

#### model 1: no pooling ###
no_pool <- "model{
    # likelihood
    for (i in 1:N){
            y[i] ~ dnorm(mu[i],tau.y)
            mu[i] <- alpha[site[i]]+bs[site[i]]*X[i,1]+bi[site[i]]*X[i,2]+
            bh[site[i]]*X[i,3]+bx[site[i]]*X[i,1]*X[i,2]
    }
    # priors for site-level parameters
    # site-level betas sampled from independent distributions
    for(j in 1:nsite){
       alpha[j]~dnorm(0,.001)
       bh[j] ~ dnorm(0,.001)
       bi[j] ~ dnorm(0,.001)
       bs[j] ~ dnorm(0,.001)
       bx[j] ~ dnorm(0,.001)
     }
    #prior for precision of y (land change rate)
    sigma.y~dunif(0,100)
    tau.y<-pow(sigma.y,-2)
 }"

#### model 2: coastwide-scale pooling ####

coast_pool <- "model{
    # likelihood
    for (i in 1:N){
            y[i] ~ dnorm(mu[i],tau.y)
            mu[i] <- alpha[site[i]]+bs[site[i]]*X[i,1]+bi[site[i]]*X[i,2]+
            bh[site[i]]*X[i,3]+bx[site[i]]*X[i,1]*X[i,2]
    }
    # priors for site-level parameters
    # site-level betas sampled from coastwide-level distribution
    for(j in 1:nsite){
       alpha[j]~dnorm(0,.001)
       bh[j] ~ dnorm(b.coast[1],tau.coast[1])
       bi[j] ~ dnorm(b.coast[2],tau.coast[2])
       bs[j] ~ dnorm(b.coast[3],tau.coast[3])
       bx[j] ~ dnorm(b.coast[4],tau.coast[4])
    }
    #prior for coastwide betas and taus
    for(j in 1:4){
      b.coast[j]~dnorm(0,.001)
      sigma.coast[j]~dunif(0,100)
      tau.coast[j]<-pow(sigma.coast[j],-2)
   }
    #prior for precision of y
    sigma.y~dunif(0,100)
    tau.y<-pow(sigma.y,-2)

 }"

#### model 3: site-basin-coastwide pooling ####
basin_coast_pool <- "model{
    # likelihood
    for (i in 1:N){
            y[i] ~ dnorm(mu[i],tau.y)
            mu[i] <- alpha[site[i]] + bs[site[i]]*X[i,1]+bi[site[i]]*X[i,2]+
            bh[site[i]]*X[i,3]+bx[site[i]]*X[i,1]*X[i,2]
    }
    # priors for site-level parameters
    # site-level betas sampled from basin-level distribution
    for(j in 1:nsite){
       alpha[j]~dnorm(0,.001)
       bh[j] ~ dmnorm(b.basin.h[basin[j]],tau.basin[1])
       bi[j] ~ dmnorm(b.basin.i[basin[j]],tau.basin[2])
       bs[j] ~ dmnorm(b.basin.s[basin[j]],tau.basin[3])
       bx[j] ~ dmnorm(b.basin.x[basin[j]],tau.basin[4])
     }
    #prior for precision of y
    tau_y~dgamma(.1,.1)
    
  #prior for site-level precision of betas
    for(i in 1:4){
      sigma.basin[i]~dunif(0,100)
      tau.basin[i]<-pow(sigma.basin[i],-2)
    }
  # basin-level betas sampled from a coastwide distribution
    for(j in 1:Nbas){
      #b.basin.h[j]~dnorm(b.coast[1],tau.coast[1])
      #b.basin.i[j]~dnorm(b.coast[2],tau.coast[2])
      #b.basin.s[j]~dnorm(b.coast[3],tau.coast[3])
      #b.basin.x[j]~dnorm(b.coast[4],tau.coast[4])
      b.basin.h[j]~dnorm(0,.001)
      b.basin.i[j]~dnorm(0,.001)
      b.basin.s[j]~dnorm(0,.001)
      b.basin.x[j]~dnorm(0,.001)
    }
  #prior for coastwide betas
   #for(j in 1:4){
   #   b.coast[j]~dnorm(0,.001)
   #   sigma.coast[j]~dunif(0,100)
   #   tau.coast[j]<-pow(sigma.coast[j],-2)
   #}
   #prior for precision of y
    sigma.y~dunif(0,100)
    tau.y<-pow(sigma.y,-2)
    
    for(x in 1:lsim){
    c.low[x]<-b.basin.s[1]*-1+b.basin.i[1]*sim[x]-1*sim[x]*b.basin.x[1]
    c.high[x]<-b.basin.s[1]*1+b.basin.i[1]*sim[x]+1*sim[x]*b.basin.x[1]
    d.low[x]<-b.basin.s[2]*-1+b.basin.i[2]*sim[x]-1*sim[x]*b.basin.x[2]
    d.high[x]<-b.basin.s[2]*1+b.basin.i[2]*sim[x]+1*sim[x]*b.basin.x[2]
    
    c.low2[x]<-b.basin.s[1]*sim[x]+b.basin.i[1]*(-1)-1*sim[x]*b.basin.x[1]
    c.high2[x]<-b.basin.s[1]*sim[x]+b.basin.i[1]*(1)+1*sim[x]*b.basin.x[1]
    d.low2[x]<-b.basin.s[2]*sim[x]+b.basin.i[2]*(-1)-1*sim[x]*b.basin.x[2]
    d.high2[x]<-b.basin.s[2]*sim[x]+b.basin.i[2]*(1)+1*sim[x]*b.basin.x[2]
    }
}"


#### set up data lists for each model####

#models 1 and 2 take the same input
dat_mod1_mod2 <- list(
  X = data.frame(df1[, c("svi", "imp", "hurr")]),
  y = df1$land_change,
  N = dim(df1)[1],
  nsite = length(unique(df1$site)),
  site = df1$site
)

#model 3 requires basin-scale info
dat_mod3 <- list(
  X = data.frame(df1[, c(5, 3, 4)]),
  y = df1$land_change,
  N = dim(df1)[1],
  nsite = length(unique(df1$site)),
  site = df1$site,
  basin = as.numeric(factor(prov_by_site)),
  Nbas = 2,
  sim=seq(-1,1,by=.1),
  lsim=length(seq(-1,1,by=.1))
)

