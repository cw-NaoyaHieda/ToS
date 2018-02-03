#局度変換を伴うsinh-arcsinh分布の単純モンテカルロ法によってVaR,ESを求める関数
SMC.fa_gairon <-function(theta){
  out <- c()
  for(i in 1:100){
    ## 正規分布からの重点サンプリングで乱数を取得
    rand.fa<-rfa_SIR_para(n=20000, mu=theta[1], 
                          sigma=theta[2], 
                          lambda=theta[3],
                          delta=theta[4])
    y <- sample(rand.fa$q,10000)
    
    # 単純モンテカルロ法
    #VaRの計算
    VaR1 <- parSapply(cl,100:10000, function(x) quantile(y[1:x], c(0.01,0.025, 0.05)))
    # ESの計算
    Es1 <- parSapply(cl,100:10000, function(x) {
      return(c( mean(y[1:x][y[1:x] < VaR1[1,x-99]]),   
                mean(y[1:x][y[1:x] < VaR1[2,x-99]]),
                mean(y[1:x][y[1:x] < VaR1[3,x-99]])
      ))})
    
    
    
    # 真値と単純モンテカルロ法の結果をまとめる
    out <- cbind(out,t(Es1))
    print(i)
  }
  return(out)
} 



# 重点サンプリングの関数(d.ISを提案分布に局度変換を伴うsinh-arcsinh分布のVaR,ESを求める)
IS.fa_gairon <- function(){
  f <- function(x, theta, par){
    exp(theta*x)*dfas2(x, mu=par[1], sigma=par[2],
                       lambda=par[3], delta=par[4])
  }
  #weightを計算するためにMを計算する．99%,97.5%,95%をまとめて行う
  M1 <- integrate(f, -30, 30, theta=theta.val1, par=fit$par2)$value
  M25 <- integrate(f, -30, 30, theta=theta.val25, par=fit$par2)$value
  M5 <- integrate(f, -30, 30, theta=theta.val5, par=fit$par2)$value
  
  #weightを計算する
  w1 <- exp(-theta.val1*rfa1)*M1
  w25 <- exp(-theta.val25*rfa25)*M25
  w5 <- exp(-theta.val5*rfa5)*M5
  
  N <-10000
  
  #99%点での計算 100~10000までサンプル数を増やして行う
  out1<- parSapply(cl, 100:N, function(x){
    # サンプルの対応するweightをくっつける
    out1<-cbind( rfa1[1:x],  w1[1:x]/x)
    # サンプルの小さい順にならべかえる
    A <- out1[sort.list(out1[,1]),]
    # weightの累積和を並べる
    A <- cbind(A, cumsum(A[,2]))
    # 累積和が0.01に一番近いサンプルが99%VaR
    # v1までのサンプルからES0.01の推定値を求める
    es1<- sum(apply(A[1:which.min(abs(A[,3]-0.01)),1:2],1,prod))/0.01
    return(es1)})
  
  #matplot(t(out1),type="l")
  #abline(h=VaR.true.FA[1],col=1)
  #abline(h=ES.true.FA[1],col=2)
  
  out25<-parSapply(cl,100:N, function(x){
    out1<-cbind(rfa25[1:x],  w25[1:x]/x)
    A <- out1[sort.list(out1[,1]),]
    A <- cbind(A, cumsum(A[,2]))
    es1<- sum(apply(A[1:which.min(abs(A[,3]-0.025)),1:2],1,prod))/0.025
    return(es1)})
  
  #matplot(t(out25),type="l")
  #abline(h=VaR.true.FA[2],col=1)
  #abline(h=ES.true.FA[2],col=2)
  
  out5<-parSapply(cl,100:N, function(x){
    out1<-cbind( rfa5[1:x],  w5[1:x]/x)
    A <- out1[sort.list(out1[,1]),]
    A <- cbind(A, cumsum(A[,2]))
    es1<- sum(apply(A[1:which.min(abs(A[,3]-0.05)),1:2],1,prod))/0.05
    return(es1)})
  
  #matplot(t(out5),type="l")
  #abline(h=VaR.true.FA[3],col=1)
  #abline(h=ES.true.FA[3],col=2)
  
  
  return(out = cbind(out1,out25,out5))
  
  }
