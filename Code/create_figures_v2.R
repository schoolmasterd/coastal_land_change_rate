
#This code is used to created figures from the samples of model 2
#load rjags library to read the saved coda-class posterior distributions
library(rjags)

#load site-level estimates
beta_samps_site<-readRDS("Output/Results/mod3_site_samples_beta.rds")
alpha_samps_site<-readRDS("Output/Results/mod3_site_samples_intercepts.rds")
prov_by_site<-readRDS("Output/Results/prov_by_site.RData")

#load province-based estimates
samps_b_prov<-readRDS("Output/Results/mod3_prov_samples.rds")
alpha_summary<-summary(alpha_samps_site)
alpha_prov<-tapply(alpha_summary$quantiles[,"50%"],prov_by_site,mean)

#load data from analysis
df1 <- read.csv("Data/combined_data_for_analysis.csv")

mod3_preds<-readRDS("Output/Results/mod3_mu_samples.rds")
#load coastal-level estimates
#samps_b_coast<-readRDS("Output/Results/mod2_coast_samples.rds")

#load site and basin names
bas_by_site<-readRDS("Output/Results/basin_by_site.rds")
site_list<-readRDS("Output/Results/site_list.rds")
basin_list<-readRDS("Output/Results/basin_list.rds")

b_prov <- summary(samps_b_prov)
1-pnorm(0,b_prov$statistics[7,1],b_prov$statistics[7,2])
1-pnorm(0,b_prov$statistics[8,1],b_prov$statistics[8,2])
####Figure 4 ####

pdf("Output/Fig_4.pdf")
par(mgp=c(3,2,0))
plot(
  c(.9,1.9,2.9,3.9),
    b_prov$quantiles[c(1,3,5,7), "50%"],
  xaxt = 'n',
  ylim = c(-1, 2),
  bty = "l",
  ylab = "Geomorphic Province-Level Mean (ha/yr)",
  xlab = "Parameter",
  xlim = c(-.2, 5),
  cex.axis=1
)
arrows(c(.9,1.9,2.9,3.9), b_prov$quantiles[c(1,3,5,7), "50%"], c(.9,1.9,2.9,3.9), b_prov$quantiles[c(1,3,5,7), "2.5%"], length = 0)
arrows(c(.9,1.9,2.9,3.9), b_prov$quantiles[c(1,3,5,7),"50%"], c(.9,1.9,2.9,3.9), b_prov$quantiles[c(1,3,5,7), "97.5%"], length = 0)
points(c(.9,1.9,2.9,3.9),
       b_prov$quantiles[c(1,3,5,7), "50%"],
       pch = 21,
       bg = "#CCCCCC")
points(c(1.1,2.1,3.1,4.1),
       b_prov$quantiles[c(2,4,6,8), "50%"],
       pch = 21,
       bg = "#333333")
arrows(c(1.1,2.1,3.1,4.1), b_prov$quantiles[c(2,4,6,8), "50%"], c(1.1,2.1,3.1,4.1), b_prov$quantiles[c(2,4,6,8), "2.5%"], length = 0)
arrows(c(1.1,2.1,3.1,4.1), b_prov$quantiles[c(2,4,6,8),"50%"], c(1.1,2.1,3.1,4.1), b_prov$quantiles[c(2,4,6,8), "97.5%"], length = 0)
abline(h = 0, lty = 2)
axis(
  side = 1,
  at = c(.9,1.9,2.9,3.9),
  labels = c("Hurricane", "Storminess", "PTD", "Interaction")
)
legend("topright",legend = c("Chenier Plain","Delta Plain"),pch=21,
       pt.bg=c("#CCCCCC","#333333"),bty='n',horiz = F,cex=1.2)
dev.off()

####Figure 5: Interaction plot####
ixn_plot<-readRDS("Output/Results/mod3_interaction_plot.rds")
inx_sum<-summary(ixn_plot)

pdf("Output/Fig_5.pdf",height = 10,width = 5)
par(mfrow=c(2,1))

ptd_raw<-(0.4719979+seq(-1,1,by=.5)*0.3023007)
plot(1:21,alpha_prov[1]+inx_sum$quantiles[grep("c.low",rownames(inx_sum$quantiles)),"50%"],type = "l",
     bty="l",xaxt="n",ylab="Land Change Rate (ha/yr)",xlab="Percent Time Drained",main="Chenier Plain",ylim=c(-2,3.5))

polygon(c(1:21,rev(1:21)),alpha_prov[1]+c(inx_sum$quantiles[grep("c.low2",rownames(inx_sum$quantiles)),"2.5%"],
                           rev(inx_sum$quantiles[grep("c.low2",rownames(inx_sum$quantiles)),"97.5%"])),
        col = adjustcolor("lightgrey",alpha.f=0.5),border=NA)
lines(1:21,alpha_prov[1]+inx_sum$quantiles[grep("c.low2",rownames(inx_sum$quantiles)),"50%"],type = "l",lwd=2)
axis(side = 1,at = c(1,6,11,16,21),labels = round(ptd_raw,1))

polygon(c(1:21,rev(1:21)),alpha_prov[1]+c(inx_sum$quantiles[grep("c.high2",rownames(inx_sum$quantiles)),"2.5%"],
                                          rev(inx_sum$quantiles[grep("c.high2",rownames(inx_sum$quantiles)),"97.5%"])),
        col=adjustcolor("lightblue",alpha.f=0.5),border = NA)
lines(1:21,alpha_prov[1]+inx_sum$quantiles[grep("c.high2",rownames(inx_sum$quantiles)),"50%"],type = "l",lty=2,lwd=2)
#Delta
mtext("(a)",side = 3,adj=0,padj=-1)
plot(1:21,alpha_prov[2]+inx_sum$quantiles[grep("d.low2",rownames(inx_sum$quantiles)),"50%"],type = "l",bty="l",xaxt="n",
     ylab="",xlab="Percent Time Drained",main = "Delta Plain",ylim=c(-2,3.5))
polygon(c(1:21,rev(1:21)),alpha_prov[2]+c(inx_sum$quantiles[grep("d.low2",rownames(inx_sum$quantiles)),"2.5%"],
                                          rev(inx_sum$quantiles[grep("d.low2",rownames(inx_sum$quantiles)),"97.5%"])),
        col = adjustcolor("lightgrey",alpha.f=0.5),border=NA)
lines(1:21,alpha_prov[2]+inx_sum$quantiles[grep("d.low",rownames(inx_sum$quantiles)),"50%"],type = "l",lwd=2)

polygon(c(1:21,rev(1:21)),alpha_prov[2]+c(inx_sum$quantiles[grep("d.high2",rownames(inx_sum$quantiles)),"2.5%"],
                                          rev(inx_sum$quantiles[grep("d.high2",rownames(inx_sum$quantiles)),"97.5%"])),
        col=adjustcolor("lightblue",alpha.f=0.5),border = NA)
lines(1:21,alpha_prov[2]+inx_sum$quantiles[grep("d.high",rownames(inx_sum$quantiles)),"50%"],type = "l",lty=2,lwd=2)
axis(side = 1,at = c(1,6,11,16,21),labels = round(ptd_raw,1))
mtext("(b)",side = 3,adj=0,padj=-1)
legend("topright",c("Calm (8.23 pa)","Stormy (149.44 pa)"),title="Storminess",lty=c(1,2),bty="n",lwd=2,xpd = T)
dev.off()

####create plots site-level predictions for appendix####
zz_heir<-summary(mod3_preds)
pretty_basin_names<-c("PO"="Pontchartrain","BS"="Breton Sound","MR"="Mississippi River Delta",
                      "BA"="Barataria","TE"="Terrebonne","AT"="Atchafalaya Delta",
                      "TV"="Teche-Vermilion","ME"="Mermentau","CS"="Calcasieu-Sabine")

cor.fit<-rep(0,length(site_list))
for(i in 1:length(site_list)){
  #set some values for plotting
  site_name<-unique(df1$site[which(df1$site==site_list[i])])
  main_title<-paste0(pretty_basin_names[names(bas_by_site[which(site_list==site_name)])]," : ",site_name)
  vals<-c(df1$land_change[which(df1$site==site_list[i])],
          zz_heir$quantiles[which(df1$site==site_list[i]),"2.5%"],zz_heir$quantiles[which(df1$site==site_list[i]),"97.5%"])
  ylims<-c(min(vals)-.5,max(vals)+.5)
  #create the plots
  pdf(paste0("Output/mod3_site_predictions_",site_name,".pdf"))
  
  plot(df1$year[which(df1$site==site_list[i])],df1$land_change[which(df1$site==site_list[i])],
       ylab="land change rate (ha/yr)",xlab="year",main=main_title,bty="l",pch=21,bg="grey",ylim=ylims)
  lines(df1$year[which(df1$site==site_list[i])],zz_heir$quantiles[which(df1$site==site_list[i]),"50%"],type="l",pch=20)
  polygon(c(df1$year[which(df1$site==site_list[i])],rev(df1$year[which(df1$site==site_list[i])])),
          c(zz_heir$quantiles[which(df1$site==site_list[i]),"2.5%"],rev(zz_heir$quantiles[which(df1$site==site_list[i]),"97.5%"])),lty=2)
  cor.fit[i]<-cor(df1$land_change[which(df1$site==site_list[i])],zz_heir$quantiles[which(df1$site==site_list[i]),"50%"])

  legend("topleft",
         legend=paste("cor = ",round(cor.fit[i],3)),bty='n')
  dev.off() 
}


###Figure 3: goodness of fit###
pdf("Output/Fig_3.pdf",height = 10,width = 6)
par(mfrow=c(2,1),oma=c(.5,.5,0,0))
###goodness of fit###
###site###
hist(cor.fit,xlab="Site-Level Correlation(Observed,Predicted)",main="")
mtext("(a)",side = 3,adj=-.15)
#cor(df1$land_change,zz_heir$quantiles[,"50%"])^2
plot(zz_heir$quantiles[,"50%"],df1$land_change,ylab="Observed Land Change Rate (ha/yr)",
     xlab="Predicted Land Change Rate (ha/yr)",pch=21,bty="l",bg=c("grey","black")[as.factor(prov_by_site)])
r2=round(1-(var(df1$land_change-zz_heir$quantiles[,"50%"]))/var(df1$land_change),3)
legend("topleft",legend = c("Chenier Plain","Delta Plain"),pch=21,
       pt.bg=c("#CCCCCC","#333333"),bty='n',horiz = F,cex=1)
mtext("(b)",side=3,adj=-.15)
text(-3,20,bquote(R^2 ==.(r2)))
dev.off()

####site estimate by basin####
b_sum<-summary(beta_samps_site)
h_who <- grep("bh", rownames(b_sum$quantiles))
i_who <- grep("bi", rownames(b_sum$quantiles))
s_who <- grep("bs", rownames(b_sum$quantiles))
x_who <- grep("bx", rownames(b_sum$quantiles))

par_group<-c("h_who","i_who","s_who","x_who")
y_labs<-c("Hurricane Estimate","Storminess Estimate","Percent Time Drained Estimate","Interaction Estimate")
pretty_basin_names<-c("PO"="Pontchartrain","BS"="Breton Sound","MR"="Mississippi River Delta",
                      "BA"="Barataria","TE"="Terrebonne","AT"="Atchafalaya Delta",
                      "TV"="Teche-Vermilion","ME"="Mermentau","CS"="Calcasieu-Sabine")



for(i in 1:9){
  pdf(paste0("Output/estimate_by_site_bas_",names(pretty_basin_names)[i],"_2025.pdf"),
      height = 10,
      width = 5)
  
  par(
    mar = c(1, 4, 2, 2),
    mfrow = c(4, 1),
    oma = c(7, 4, 3, 1)
  )
  for(j in 1:4){
    par_to_plot<-get(par_group[j])
    ylim_<-range(c(b_sum$quantiles[par_to_plot[bas_by_site==i], "50%"],
                   b_sum$quantiles[par_to_plot[bas_by_site==i], "2.5%"],
                   b_sum$quantiles[par_to_plot[bas_by_site==i], "97.5%"]))*c(.95,1.05)
    
    plot(
      1:length(par_to_plot[bas_by_site==i]),
      b_sum$quantiles[par_to_plot[bas_by_site==i], "50%"],
      ylim = ylim_,
      pch = 21,
      bg = "darkgrey",
      ylab = y_labs[j],
      bty = "l",
      xaxt = "n",
      xlab = ""
    )
    arrows(1:length(par_to_plot[bas_by_site==i]),
           b_sum$quantiles[par_to_plot[bas_by_site==i], "50%"],
           1:length(par_to_plot[bas_by_site==i]),
           b_sum$quantiles[par_to_plot[bas_by_site==i], "2.5%"],
           length = 0)
    arrows(1:length(par_to_plot[bas_by_site==i]),
           b_sum$quantiles[par_to_plot[bas_by_site==i], "50%"],
           1:length(par_to_plot[bas_by_site==i]),
           b_sum$quantiles[par_to_plot[bas_by_site==i], "97.5%"],
           length = 0)
    
    abline(h=0,lwd=1,lty=2)
  }
  
  h_nms <-
    site_list[h_who[bas_by_site==i]]
  axis(
    side = 1,
    at = 1:length(h_who[bas_by_site==i]),
    labels = h_nms,
    las = 2,
    cex.axis =.5
  )
  mtext(pretty_basin_names[i],side = 3,outer = T)
  mtext("Site", side = 1, outer=T,padj = 6)
  dev.off()
}
####plot alphas ####
#samps_alpha_heir<-readRDS("Output/Results/site_level_intercepts.rds")
alphs<-summary(alpha_samps_site)
dim(alphs$quantiles)

plt_matrix<-matrix(c(1,2,3,4,5,6,7,8,9),nrow=3,ncol=3)
for(i in 1:3){
  pdf(paste0("Output/estimate_by_site_bas_alpha_",i,".pdf"),
      height = 10,
      width = 5)
  
  par(
    mar = c(6, 4, 2, 2),
    mfrow = c(3, 1),
    oma = c(5, 4, 3, 1)
  )
  for(j in 1:3){
    par_to_plot<-alphs$quantiles[bas_by_site==plt_matrix[j,i],]
    ylim_<-range(c(par_to_plot[, "50%"],
                   par_to_plot[, "2.5%"],
                   par_to_plot[, "97.5%"]))*c(.95,1.05)
    len<-dim(par_to_plot)[1]
    plot(
      1:len,
      par_to_plot[, "50%"],
      ylim = ylim_,
      pch = 21,
      bg = "darkgrey",
      ylab = "Intercept Estimate",
      bty = "l",
      xaxt = "n",
      xlab = "",
      main = pretty_basin_names[plt_matrix[j,i]]
    )
    arrows(1:len,
           par_to_plot[, "50%"],
           1:len,
           par_to_plot[, "2.5%"],
           length = 0)
    arrows(1:len,
           par_to_plot[, "50%"],
           1:len,
           par_to_plot[, "97.5%"],
           length = 0)
    
    abline(h=0,lwd=1,lty=2)
    
    h_nms <-
      site_list[bas_by_site==plt_matrix[j,i]]
    axis(
      side = 1,
      at = 1:len,
      labels = h_nms,
      las = 2,
      cex.axis =.5
    )
  }
  #mtext(pretty_basin_names[i],side = 3,outer = T)
  mtext("Site", side = 1, outer=T,padj = 1)
  dev.off()
}
