cl <- makeCluster(detectCores()-1)  # クラスタの作成
clusterExport(cl,c('s_inverse','dfas2','dfas','Cf','Sf'))


#mu.val1 <- ES.true.FA[1]
#mu.val25 <- ES.true.FA[2]
#mu.val5 <- ES.true.FA[3]

N <-10000
### 次の関数を新しく用意した。
rIS_SIR2_ <- function(n, par, par2){
  ## 正規分布を提案分布に par2 = c( mu, sigma) 
  q <- rnorm(n, mean=par2[1], sd=par2[2])
  ## 重み
  #hosts <- rep('localhost',12)
  #hosts <- 7
  #cl <- makeCluster(hosts, "SOCK")
  w <- parSapply(cl,q, 
                 dfas2, mu = par[1], sigma= par[2] , lambda=par[3], delta = par[4]
  )/dnorm(q, mean=par2[1], sd=par2[2])
  w <- w/sum(w)
  ## resample
  q.resample <- Resample1(q, weight=w, NofSample = n)
  list( q=q.resample, w=w)
}


out <- c()
rfa.IS.1<-rIS_SIR2_(n=2 * N, par=c( mu.val1, fit$par2[2:4]),
                    par2 = c(mu.val1, sd(rt))) 

rfa.IS.25<-rIS_SIR2_(n=2 * N, par=c( mu.val25, fit$par2[2:4]),
                     par2 = c(mu.val25, sd(rt))) 

rfa.IS.5<-rIS_SIR2_(n=2 * N, par=c( mu.val5, fit$par2[2:4]),
                    par2 = c(mu.val5, sd(rt))) 

rfa1 <- sample(rfa.IS.1$q, N)
rfa25 <- sample(rfa.IS.25$q, N)
rfa5 <- sample(rfa.IS.5$q, N)

clusterExport(cl,c('s_inverse','dfas2','dfas','Cf','Sf',
                   'rfa1','rfa25','rfa5'))  
# 重点サンプリングの関数の関数も新しいものにした。
f <- function(x, par, par2){
  #par2 f # par g
  dfas2(x, mu=par[1], sigma=par[2],
        lambda=par[3], delta=par[4])/
    dfas2(x, mu=par2[1], sigma=par[2],
          lambda=par[3], delta=par[4])}

#weightを計算する
w1 <- f(rfa1, par = fit$par2, par2=c(mu.val1, fit$par2[2:4]))
w25 <- f(rfa25, par = fit$par2, par2=c(mu.val25, fit$par2[2:4]))
w5 <- f(rfa5, par = fit$par2, par2=c(mu.val5, fit$par2[2:4]))


#hosts <- rep('localhost',12)
#hosts <- 7
#cl <- makeCluster(hosts, "SOCK")
clusterExport(cl,c('w1','w25','w5'))  
#99%点での計算 100~10000までサンプル数を増やして行う
out1<- parSapply(cl, 100:N, function(x){
  # サンプルの対応するweightをくっつける
  out1<-cbind(rfa1[1:x],  w1[1:x]/x)
  # サンプルの小さい順にならべかえる
  A <- out1[sort.list(out1[,1]),]
  # weightの累積和を並べる
  A <- cbind(A, cumsum(A[,2]))
  # 累積和が0.01に一番近いサンプルが99%VaR
  v1 <- A[which.min(abs(A[,3]-0.01)),1]
  # v1までのサンプルからES0.01の推定値を求める
  if(is.null(dim(A[1:which.min(abs(A[,3]-0.01)),]))){
    es1 <- sum(prod(A[1:which.min(abs(A[,3]-0.01)),c(1:2)]))/0.01
  }else{
    es1 <- sum(apply(A[1:which.min(abs(A[,3]-0.01)),1:2],1,prod))/0.01
  }
  return(c(v1, es1))})

out25<-parSapply(cl, 100:N, function(x){
  out1<-cbind( rfa25[1:x],  w25[1:x]/x)
  A <- out1[sort.list(out1[,1]),]
  A <- cbind(A, cumsum(A[,2]))
  v1<-A[which.min(abs(A[,3]-0.025)),1]
  if(is.null(dim(A[1:which.min(abs(A[,3]-0.025)),]))){
    es1 <- sum(prod(A[1:which.min(abs(A[,3]-0.025)),c(1:2)]))/0.025
  }else{
    es1<- sum(apply(A[1:which.min(abs(A[,3]-0.025)),1:2],1,prod))/0.025
  }
  return(c(v1, es1))})


out5<-parSapply(cl, 100:N, function(x){
  out1<-cbind( rfa5[1:x],  w5[1:x]/x)
  A <- out1[sort.list(out1[,1]),]
  A <- cbind(A, cumsum(A[,2]))
  v1<-A[which.min(abs(A[,3]-0.05)),1]
  if(is.null(dim(A[1:which.min(abs(A[,3]-0.05)),]))){
    es1 <- sum(prod(A[1:which.min(abs(A[,3]-0.05)),c(1:2)]))/0.05
  }else{
    es1<- sum(apply(A[1:which.min(abs(A[,3]-0.05)),1:2],1,prod))/0.05
  }
  return(c(v1, es1))})

out <- rbind(out1, out25, out5)

while(dim(out)[1] < 601){
  rfa.IS.1<-rIS_SIR2_(n=2 * N, par=c( mu.val1, fit$par2[2:4]),
                      par2 = c(mu.val1, sd(rt))) 
  
  rfa.IS.25<-rIS_SIR2_(n=2 * N, par=c( mu.val25, fit$par2[2:4]),
                       par2 = c(mu.val25, sd(rt))) 
  
  rfa.IS.5<-rIS_SIR2_(n=2 * N, par=c( mu.val5, fit$par2[2:4]),
                      par2 = c(mu.val5, sd(rt))) 
  
  rfa1 <- sample(rfa.IS.1$q, N)
  rfa25 <- sample(rfa.IS.25$q, N)
  rfa5 <- sample(rfa.IS.5$q, N)
  
  clusterExport(cl,c('s_inverse','dfas2','dfas','Cf','Sf',
                     'rfa1','rfa25','rfa5'))  
  # 重点サンプリングの関数の関数も新しいものにした。
  f <- function(x, par, par2){
    #par2 f # par g
    dfas2(x, mu=par[1], sigma=par[2],
          lambda=par[3], delta=par[4])/
      dfas2(x, mu=par2[1], sigma=par[2],
            lambda=par[3], delta=par[4])}
  
  #weightを計算する
  w1 <- f(rfa1, par = fit$par2, par2=c(mu.val1, fit$par2[2:4]))
  w25 <- f(rfa25, par = fit$par2, par2=c(mu.val25, fit$par2[2:4]))
  w5 <- f(rfa5, par = fit$par2, par2=c(mu.val5, fit$par2[2:4]))
  

  #hosts <- rep('localhost',12)
  #hosts <- 7
  #cl <- makeCluster(hosts, "SOCK")
  clusterExport(cl,c('w1','w25','w5'))  
  #99%点での計算 100~10000までサンプル数を増やして行う
  out1<- parSapply(cl, 100:N, function(x){
    # サンプルの対応するweightをくっつける
    out1<-cbind(rfa1[1:x],  w1[1:x]/x)
    # サンプルの小さい順にならべかえる
    A <- out1[sort.list(out1[,1]),]
    # weightの累積和を並べる
    A <- cbind(A, cumsum(A[,2]))
    # 累積和が0.01に一番近いサンプルが99%VaR
    v1 <- A[which.min(abs(A[,3]-0.01)),1]
    # v1までのサンプルからES0.01の推定値を求める
    if(is.null(dim(A[1:which.min(abs(A[,3]-0.01)),]))){
      es1 <- sum(prod(A[1:which.min(abs(A[,3]-0.01)),c(1:2)]))/0.01
    }else{
      es1 <- sum(apply(A[1:which.min(abs(A[,3]-0.01)),1:2],1,prod))/0.01
    }
    return(c(v1, es1))})
  
  out25<-parSapply(cl, 100:N, function(x){
    out1<-cbind( rfa25[1:x],  w25[1:x]/x)
    A <- out1[sort.list(out1[,1]),]
    A <- cbind(A, cumsum(A[,2]))
    v1<-A[which.min(abs(A[,3]-0.025)),1]
    if(is.null(dim(A[1:which.min(abs(A[,3]-0.025)),]))){
      es1 <- sum(prod(A[1:which.min(abs(A[,3]-0.025)),c(1:2)]))/0.025
    }else{
      es1<- sum(apply(A[1:which.min(abs(A[,3]-0.025)),1:2],1,prod))/0.025
    }
    return(c(v1, es1))})
  
  
  out5<-parSapply(cl, 100:N, function(x){
    out1<-cbind( rfa5[1:x],  w5[1:x]/x)
    A <- out1[sort.list(out1[,1]),]
    A <- cbind(A, cumsum(A[,2]))
    v1<-A[which.min(abs(A[,3]-0.05)),1]
    if(is.null(dim(A[1:which.min(abs(A[,3]-0.05)),]))){
      es1 <- sum(prod(A[1:which.min(abs(A[,3]-0.05)),c(1:2)]))/0.05
    }else{
      es1<- sum(apply(A[1:which.min(abs(A[,3]-0.05)),1:2],1,prod))/0.05
    }
    return(c(v1, es1))})
  
  out <- rbind(out, out1, out25, out5)
}
stopCluster(cl)


out <- t(out)
IS_VaR_1 <- out[,c(seq(1,600,6))]
IS_ES_1 <- out[,c(seq(2,600,6))]
IS_VaR_25 <- out[,c(seq(3,600,6))]
IS_ES_25 <- out[,c(seq(4,600,6))]
IS_VaR_5 <- out[,c(seq(5,600,6))]
IS_ES_5 <- out[,c(seq(6,600,6))]

IS_VaR_mean <- data.frame(rowMeans(IS_VaR_1),rowMeans(IS_VaR_25),rowMeans(IS_VaR_5))
IS_VaR_sd <- data.frame(apply(IS_VaR_1,1,sd),apply(IS_VaR_25,1,sd),apply(IS_VaR_5,1,sd))

IS_ES_mean <- data.frame(rowMeans(IS_ES_1),rowMeans(IS_ES_25),rowMeans(IS_ES_5))
IS_ES_sd <- data.frame(apply(IS_ES_1,1,sd),apply(IS_ES_25,1,sd),apply(IS_ES_5,1,sd))


range_95_VaR_IS_min <- data.frame(100:N, -IS_VaR_mean + 2 * IS_VaR_sd)
range_95_VaR_IS_max <- data.frame(100:N, -IS_VaR_mean - 2 * IS_VaR_sd)
range_95_ES_IS_min <- data.frame(100:N, -IS_ES_mean + 2 * IS_ES_sd)
range_95_ES_IS_max <- data.frame(100:N, -IS_ES_mean - 2 * IS_ES_sd)

names(range_95_VaR_IS_min)<- c("N","VaR99","VaR97.5","VaR95")
names(range_95_VaR_IS_max)<- c("N","VaR99","VaR97.5","VaR95")
names(range_95_ES_IS_min)<- c("N","ES99","ES97.5","ES95")
names(range_95_ES_IS_max)<- c("N","ES99","ES97.5","ES95")
range_95_VaR_IS_min=melt(range_95_VaR_IS_min,id.vars=c("N"))
range_95_VaR_IS_max=melt(range_95_VaR_IS_max,id.vars=c("N"))
range_95_ES_IS_min=melt(range_95_ES_IS_min,id.vars=c("N"))
range_95_ES_IS_max=melt(range_95_ES_IS_max,id.vars=c("N"))


range_95_min <- rbind(range_95_VaR_IS_min,range_95_ES_IS_min)
range_95_max <- rbind(range_95_VaR_IS_max,range_95_ES_IS_max)

out1 <-data.frame(100:N, -IS_VaR_mean)
out2 <-data.frame(100:N, -IS_ES_mean)
names(out1)<- c("N","VaR99","VaR97.5","VaR95")
names(out2)<- c("N","ES99","ES97.5","ES95")

df1=melt(out1,id.vars=c("N"))
df2=melt(out2,id.vars=c("N"))

df <- rbind(df1,df2)


df <- data.frame(df, max=range_95_min[, 3], min = range_95_max[, 3])

df$variable <- factor(df$variable, levels=c("VaR99","ES99","VaR97.5","ES97.5","VaR95","ES95"))

