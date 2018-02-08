#pattern3

backcolour <- c("#FDE3E1","#F0EBCC","#CCF1D7","#CCF2F3","#DFEBFF","#FDE0F9")

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
  xlim(c(0,10000))+scale_fill_manual(values = c( backcolour))

print(plot_b)


#pattern5
plot_b=ggplot(data=df)+
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
  geom_line(aes(x=N,y=min,color=variable))+
  geom_line(aes(x=N,y=max,color=variable))+
  xlim(c(0,10000))

print(plot_b)


#
plot_b=ggplot(data=df)+
  geom_ribbon(aes(x=N,ymin=min,ymax=max,fill=variable))+
  geom_hline(aes(yintercept=-ES.true.FA[1]))+
  geom_hline(aes(yintercept=-ES.true.FA[2]))+
  geom_hline(aes(yintercept=-ES.true.FA[3]))+
  geom_line(aes(x=N,y=min,color=variable))+
  geom_line(aes(x=N,y=max,color=variable))+
  geom_line(aes(N,value,lty=variable,col=variable),size=.8)+
  theme_bw()+ylab("Expected Shortfall")+
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black"),
    legend.title=element_blank(),legend.position=c(0.25,0.82),
    legend.direction = "horizontal",legend.key.width = unit(12, "points"))+
  facet_grid(.~method)+
  xlim(c(0,10000))+scale_fill_manual(values = c( backcolour))

print(plot_b)
