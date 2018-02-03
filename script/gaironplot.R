#面倒なので流用しています
library(parallel)
library(ggplot2)
source("script/functions.R")
source("script/functions_gaironplot.R")
fit<-c()
fit$par2 <- c(0.144,0.484,-0.338,0.563)

#真値計算
#99%,97.5%,95%の各点に対して，先ほどの関数を用いて求める
VaR1.fa <- qfas(0.01, mu=fit$par2[1], sigma=fit$par2[2],
                lambda=fit$par2[3], delta = fit$par2[4])
VaR25.fa <- qfas(0.025, mu=fit$par2[1], sigma=fit$par2[2],
                 lambda=fit$par2[3], delta = fit$par2[4])
VaR5.fa <- qfas(0.05, mu=fit$par2[1], sigma=fit$par2[2],
                lambda=fit$par2[3], delta = fit$par2[4])
#まとめておく
VaR.true.FA <- c(VaR1.fa ,VaR25.fa ,VaR5.fa )
ES1.fa <- find.ES(p=0.01, par=fit$par2)
ES25.fa <- find.ES(p=0.025, par=fit$par2)
ES5.fa <- find.ES(p=0.05, par=fit$par2)
#まとめる
ES.true.FA <- c(ES1.fa ,ES25.fa ,ES5.fa )

cl <- makeCluster(25, type="SOCK")
cl_l <- length(cl)
clusterExport(cl,list("s_inverse","dfas","dfas2","Cf","Sf","cl_l","fit","rfa_SIR","Resample1"))



SMC.fa.out <- SMC.fa_gairon(fit$par2)

SMC_ES_1 <- SMC.fa.out[,c(seq(1,300,3))]
SMC_ES_25 <- SMC.fa.out[,c(seq(2,300,3))]
SMC_ES_5 <- SMC.fa.out[,c(seq(3,300,3))]

SMC_ES_mean <- data.frame(rowMeans(SMC_ES_1),rowMeans(SMC_ES_25),rowMeans(SMC_ES_5))
SMC_ES_sd <- data.frame(apply(SMC_ES_1,1,sd),apply(SMC_ES_25,1,sd),apply(SMC_ES_5,1,sd))


## IS(重点サンプリング)による<br>VaRとESを計算
# 99%,97.5%,95%それぞれのVaRと平均が一致するthetaを取得
theta.val1<- find.theta(0.01, fit$par2)
theta.val25<- find.theta(0.025, fit$par2)
theta.val5<- find.theta(0.05, fit$par2)


clusterExport(cl,list("VaR.true.FA"))

out.fa<-c()
rfa.IS.1<-rIS_SIR(n=20000, par=fit$par2, par2=VaR.true.FA[1], theta=theta.val1)
rfa.IS.25<-rIS_SIR(n=20000, par=fit$par2, par2=VaR.true.FA[2], theta=theta.val25)
rfa.IS.5<-rIS_SIR(n=20000, par=fit$par2, par2=VaR.true.FA[3], theta=theta.val5)
# サンプリングしたものを入力としてFA分布の重点サンプリングを行う
rfa1 <- sample(rfa.IS.1$q, 10000)
rfa25 <- sample(rfa.IS.25$q, 10000)
rfa5 <- sample(rfa.IS.5$q, 10000)
clusterExport(cl,list("rfa1","rfa25","rfa5"))
try(out.fa <- cbind(out.fa,IS.fa_gairon()))


while(dim(out.fa)[2] < 301){
  rfa.IS.1<-rIS_SIR(n=20000, par=fit$par2, par2=VaR.true.FA[1], theta=theta.val1)
  rfa.IS.25<-rIS_SIR(n=20000, par=fit$par2, par2=VaR.true.FA[2], theta=theta.val25)
  rfa.IS.5<-rIS_SIR(n=20000, par=fit$par2, par2=VaR.true.FA[3], theta=theta.val5)
  # サンプリングしたものを入力としてFA分布の重点サンプリングを行う
  rfa1 <- sample(rfa.IS.1$q, 10000)
  rfa25 <- sample(rfa.IS.25$q, 10000)
  rfa5 <- sample(rfa.IS.5$q, 10000)
  clusterExport(cl,list("rfa1","rfa25","rfa5"))
  try(out.fa <- cbind(out.fa,IS.fa_gairon()))
  print(dim(out.fa)[2])
}

stopCluster(cl)

IS_ES_1 <- out.fa[,c(seq(1,300,3))]
IS_ES_25 <- out.fa[,c(seq(2,300,3))]
IS_ES_5 <- out.fa[,c(seq(3,300,3))]

IS_ES_mean <- data.frame(rowMeans(IS_ES_1),rowMeans(IS_ES_25),rowMeans(IS_ES_5))
IS_ES_sd <- data.frame(apply(IS_ES_1,1,sd),apply(IS_ES_25,1,sd),apply(IS_ES_5,1,sd))



range_95_IS_min <- data.frame(100:10000, -IS_ES_mean + 2 * IS_ES_sd)
range_95_IS_max <- data.frame(100:10000, -IS_ES_mean - 2 * IS_ES_sd)
range_95_SMC_min <- data.frame(100:10000, -SMC_ES_mean + 2 * SMC_ES_sd)
range_95_SMC_max <- data.frame(100:10000, -SMC_ES_mean - 2 * SMC_ES_sd)
names(range_95_IS_min)<- c("N","ES99(IS)","ES97.5(IS)","ES95(IS)")
names(range_95_IS_max)<- c("N","ES99(IS)","ES97.5(IS)","ES95(IS)")
names(range_95_SMC_min)<- c("N","ES99(SMC)","ES97.5(SMC)","ES95(SMC)")
names(range_95_SMC_max)<- c("N","ES99(SMC)","ES97.5(SMC)","ES95(SMC)")
range_95_IS_min=melt(range_95_IS_min,id.vars=c("N"))
range_95_IS_max=melt(range_95_IS_max,id.vars=c("N"))
range_95_SMC_min=melt(range_95_SMC_min,id.vars=c("N"))
range_95_SMC_max=melt(range_95_SMC_max,id.vars=c("N"))

tmp1 <- cbind(range_95_IS_min, method="IS")
tmp2 <- cbind(range_95_SMC_min, method="SMC")
range_95_min <- rbind(tmp1,tmp2)

tmp1 <- cbind(range_95_IS_max, method="IS")
tmp2 <- cbind(range_95_SMC_max, method="SMC")
range_95_max <- rbind(tmp1,tmp2)

out1 <-data.frame(100:10000, -IS_ES_mean)
out2 <-data.frame(100:10000, -SMC_ES_mean)
names(out1)<- c("N","ES99(IS)","ES97.5(IS)","ES95(IS)")
names(out2)<- c("N","ES99(SMC)","ES97.5(SMC)","ES95(SMC)")

df1=melt(out1,id.vars=c("N"))
df2=melt(out2,id.vars=c("N"))

df1 <- cbind(df1, method="IS")
df2 <- cbind(df2, method="SMC")
df <- rbind(df1,df2)


df <- data.frame(df, max=range_95_min[, 3], min = range_95_max[, 3])

df$variable <- factor(df$variable, levels=c("ES99(IS)","ES99(SMC)","ES97.5(IS)","ES97.5(SMC)","ES95(IS)","ES95(SMC)"))

#pattern1
plot_b=ggplot(data=df)+
  geom_ribbon(aes(x=N,ymin=min,ymax=max,fill=variable))+
  geom_hline(aes(yintercept=-ES.true.FA[1]))+
  geom_hline(aes(yintercept=-ES.true.FA[2]))+
  geom_hline(aes(yintercept=-ES.true.FA[3]))+
  geom_line(aes(N,value,lty=variable,col=variable),size=.8)+
  theme_bw()+ylab("Expected Shortfall")+
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black"),
    legend.title=element_blank(),legend.position=c(0.25,0.82),
    legend.direction = "horizontal",legend.key.width = unit(12, "points"))+
  facet_grid(.~method)+
  xlim(c(0,10000))

print(plot_b)

#pattern2
plot_b=ggplot(data=df)+
  geom_line(aes(N,value,lty=variable,col=variable),size=.8)+
  geom_hline(aes(yintercept=-ES.true.FA[1]))+
  geom_hline(aes(yintercept=-ES.true.FA[2]))+
  geom_hline(aes(yintercept=-ES.true.FA[3]))+
  theme_bw()+ylab("Expected Shortfall")+
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black"),
    legend.title=element_blank(),legend.position=c(0.25,0.92),
    legend.direction = "horizontal",legend.key.width = unit(12, "points"))+
  facet_grid(method~.)+
  geom_ribbon(aes(x=N,ymin=min,ymax=max,fill=variable),alpha=0.2)+
  xlim(c(0,10000))


print(plot_b)

#pattern3
plot_b=ggplot(data=df)+
  geom_ribbon(aes(x=N,ymin=min,ymax=max,fill=variable),alpha=0.2)+
  geom_hline(aes(yintercept=-ES.true.FA[1]))+
  geom_hline(aes(yintercept=-ES.true.FA[2]))+
  geom_hline(aes(yintercept=-ES.true.FA[3]))+
  geom_line(aes(N,value,lty=variable,col=variable),size=.8)+
  theme_bw()+ylab("Expected Shortfall")+
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black"),
    legend.title=element_blank(),legend.position=c(0.35,0.9),
    legend.direction = "horizontal",legend.key.width = unit(12, "points"))+
  xlim(c(0,10000))






SMC_ES_1 <- SMC.fa.out[,c(seq(1,300,3))]
SMC_ES_25 <- SMC.fa.out[,c(seq(2,300,3))]
SMC_ES_5 <- SMC.fa.out[,c(seq(3,300,3))]

IS_ES_1 <- out.fa[,c(seq(1,300,3))]
IS_ES_25 <- out.fa[,c(seq(2,300,3))]
IS_ES_5 <- out.fa[,c(seq(3,300,3))]

IS_ES_1 <- data.frame(IS_ES_1 %>% melt(),name="ES99(IS)")
IS_ES_25 <- data.frame(IS_ES_25 %>% melt(),name="ES97.5(IS)")
IS_ES_5 <- data.frame(IS_ES_5 %>% melt(),name="ES95(IS)")
SMC_ES_1 <- data.frame(SMC_ES_1 %>% melt(),name="ES99(SMC)")
SMC_ES_25 <- data.frame(SMC_ES_25 %>% melt(),name="ES97.5(SMC)")
SMC_ES_5 <- data.frame(SMC_ES_5 %>% melt(),name="ES95(SMC)")




ES_box_IS <- rbind(IS_ES_1[,-2], IS_ES_25[,-2], IS_ES_5[,-2])
ES_box_SMC <- rbind(SMC_ES_1[,-2], SMC_ES_25[,-2], SMC_ES_5[,-2])
ES_box <- rbind(data.frame(ES_box_IS, method="IS"),
                data.frame(ES_box_SMC, method="SMC"))
ES_box$Var1 <- ES_box$Var1 + 99


#pattern4
ggplot(ES_box %>% filter(Var1%%500 == 0), aes(x = as.factor(Var1),y = -value,fill=name,color=name))+
  geom_boxplot(position = "identity")+geom_hline(aes(yintercept=-ES.true.FA[1]))+
  geom_hline(aes(yintercept=-ES.true.FA[2]))+
  geom_hline(aes(yintercept=-ES.true.FA[3]))+
  theme_bw()+ylab("Expected Shortfall")+
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black"))+
  xlab("N")+
  facet_grid(method~.,scales = "free")+ylim(2.5,6.2)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


  
#pattern5
plot_b=ggplot(data=df)+
  geom_line(aes(N,value,lty=variable,col=variable),size=.8)+
  geom_hline(aes(yintercept=-ES.true.FA[1]))+
  geom_hline(aes(yintercept=-ES.true.FA[2]))+
  geom_hline(aes(yintercept=-ES.true.FA[3]))+
  theme_bw()+ylab("Expected Shortfall")+
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black"),
    legend.title=element_blank(),legend.position=c(0.25,0.82),
    legend.direction = "horizontal",legend.key.width = unit(12, "points"))+
  facet_grid(.~method)+
  geom_line(aes(x=N,y=min,color=variable))+
  geom_line(aes(x=N,y=max,color=variable))+
  xlim(c(0,10000))

print(plot_b)


#pattern6
plot_b=ggplot(data=df)+
  geom_line(aes(N,value,lty=variable,col=variable),size=.8)+
  geom_hline(aes(yintercept=-ES.true.FA[1]))+
  geom_hline(aes(yintercept=-ES.true.FA[2]))+
  geom_hline(aes(yintercept=-ES.true.FA[3]))+
  theme_bw()+ylab("Expected Shortfall")+
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black"),
    legend.title=element_blank(),legend.position=c(0.35,0.9),
    legend.direction = "horizontal",legend.key.width = unit(12, "points"))+
  facet_grid(method~.)+
  geom_line(aes(x=N,y=min,color=variable))+
  geom_line(aes(x=N,y=max,color=variable))+
  xlim(c(0,10000))

print(plot_b)


#patter7
plot_b=ggplot(data=df)+
  geom_line(aes(N,value,lty=variable,col=variable),size=.8)+
  geom_hline(aes(yintercept=-ES.true.FA[1]))+
  geom_hline(aes(yintercept=-ES.true.FA[2]))+
  geom_hline(aes(yintercept=-ES.true.FA[3]))+
  theme_bw()+ylab("Expected Shortfall")+
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black"),
    legend.title=element_blank(),legend.position=c(0.35,0.9),
    legend.direction = "horizontal",legend.key.width = unit(12, "points"))+
  geom_line(aes(x=N,y=min,color=variable))+
  geom_line(aes(x=N,y=max,color=variable))+
  xlim(c(0,10000))

print(plot_b)


#pattern8
ggplot(data=df)+
  geom_errorbar(aes(x=N, ymax = max, ymin = min, color=variable),width = 0.1)+
  geom_hline(aes(yintercept=-ES.true.FA[1]))+
  geom_hline(aes(yintercept=-ES.true.FA[2]))+
  geom_hline(aes(yintercept=-ES.true.FA[3]))+
  geom_point(aes(N,value,lty=variable,col=variable),size=.8)+
  theme_bw()+ylab("Expected shortfall")+
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black"),
    legend.title=element_blank(),legend.position=c(0.32,0.9),
    legend.direction = "horizontal",legend.key.width = unit(12, "points"))+
  facet_grid(method~.)+
  xlim(c(0,10000))


#pattern9
ggplot(data=df)+
  geom_hline(aes(yintercept=-ES.true.FA[1]))+
  geom_hline(aes(yintercept=-ES.true.FA[2]))+
  geom_hline(aes(yintercept=-ES.true.FA[3]))+
  geom_line(aes(N,value,lty=variable,col=variable),size=.8)+
  theme_bw()+ylab("Expected shortfall")+
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black"),
    legend.title=element_blank(),legend.position=c(0.25,0.9),
    legend.direction = "horizontal",legend.key.width = unit(12, "points"))+
  geom_errorbar(aes(x=N, ymax = max, ymin = min, color=variable),alpha=0.01,width = 0.1)+
  facet_grid(.~method)+
  xlim(c(0,10000))

#pattern10
ggplot(data=df)+
  geom_hline(aes(yintercept=-ES.true.FA[1]))+
  geom_hline(aes(yintercept=-ES.true.FA[2]))+
  geom_hline(aes(yintercept=-ES.true.FA[3]))+
  geom_line(aes(N,value,lty=variable,col=variable),size=.8)+
  geom_point(aes(N,value,lty=variable,col=variable),size=.8)+
  theme_bw()+ylab("Expected shortfall")+
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black"),
    legend.title=element_blank(),legend.position=c(0.25,0.9),
    legend.direction = "horizontal",legend.key.width = unit(12, "points"))+
  geom_errorbar(aes(x=N, ymax = max, ymin = min, color=variable),alpha=0.01,width = 0.1)+
  facet_grid(.~method)+
  xlim(c(0,10000))

#pattern11
ggplot(data=df %>% filter(N%%50 == 0))+
  geom_hline(aes(yintercept=-ES.true.FA[1]))+
  geom_hline(aes(yintercept=-ES.true.FA[2]))+
  geom_hline(aes(yintercept=-ES.true.FA[3]))+
  theme_bw()+ylab("Expected shortfall")+
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black"),
    legend.title=element_blank(),legend.position=c(0.32,0.9),
    legend.direction = "horizontal",legend.key.width = unit(12, "points"))+
  geom_errorbar(aes(x=N, ymax = max, ymin = min, color=variable), width = 0.1)+
  geom_line(aes(N,value,lty=variable,col=variable),alpha=0.1,size=.8)+
  facet_grid(method~.)+
  xlim(c(0,10000))

#pattern12
ggplot(data=df %>% filter(N%%100 == 0))+
  geom_line(aes(N,value,lty=variable,col=variable),size=.8)+
  geom_hline(aes(yintercept=-ES.true.FA[1]))+
  geom_hline(aes(yintercept=-ES.true.FA[2]))+
  geom_hline(aes(yintercept=-ES.true.FA[3]))+
  theme_bw()+ylab("Expected shortfall")+
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black"),
    legend.title=element_blank(),legend.position=c(0.25,0.9),
    legend.direction = "horizontal",legend.key.width = unit(12, "points"))+
  geom_errorbar(aes(x=N, ymax = max, ymin = min, color=variable), width = 0.1)+
  facet_grid(.~method)+
  xlim(c(0,10000))
