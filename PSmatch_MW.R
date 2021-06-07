#1:1 PSM MW
#点推定

#install.packages("geeM_0.10.1.tar.gz", repos = NULL, type = "source")
#library、その他もろもろ
library(Matching)
library(tidyverse)
source("boot.R")
source("CI.R")

#推定量の保存
n <- 1000 #データ数
J <- 1000
Mean <- matrix(rep(0,J*2),J,2) #estimatorsに結果を代入
boot_sd <-matrix(rep(0,J*4),J,4) #ブートストラップ標準偏差
w <- rep(0,n) #wはmatching weightの重み


#logitモデル、アウトカムモデルのデータ生成
beta <- matrix(c(-1, 1, -1.2, 2, -1.5),1,5)
alpha <- matrix(c(0.5, 1, 0.6, 2.2, -1.2),1,5)
delta <- 2.5  #真の因果効果

for (j in 1:J){
  set.seed(j) #シード設定
  #データ生成
  x_0 <- rep(1, n) 
  x_4 <- rbinom(n,1,0.5) #B(0.5)に従う分布
  x_3 <- rbinom(n,1,x_4*0.6 +(1-x_4)*0.4) #B(x_4*0.6 + (1-x_4)*0.4)
  
  #x_1,x_2は多変量正規分布コレスキー分解
  #z_1,z_2の標準正規分布をもとに生成
  z_1 <- rnorm(n,0,1)
  z_2 <- rnorm(n,0,1) 
  #分散共分散行列のために必要な行列生成
  m_1 <- matrix(c(1,0.5,0.5,1),2,2)
  m_2 <- matrix(c(2,0.25,0.25,2),2,2)
  
  x_1 <- rep(0,n) 
  x_2 <- rep(0,n)
  t1 <- mean(-x_3+ x_4 + 0.5*x_3*x_4)
  t2 <- mean(x_3-x_4+x_3*x_4)
  for(i in 1:n){
    cov_x3 <- x_3[i]*m_1 + (1-x_3[i])*m_2 #共分散行列作成
    #コレスキー分解 chol(cov_x3)
    #xに変換した変数を入れている L*z: Lが下三角行列
    x <- t(chol(cov_x3))%*%c(z_1[i],z_2[i])+c(t1,t2)
      #c(-x_3[i]+x_4[i]+0.5*x_3[i]*x_4[i],x_3[i]-x_4[i]+x_3[i]*x_4[i])
    x_1[i] <- x[1]
    x_2[i] <- x[2]
  }
  
  #x_0~x_4をまとめる(Data)
  Data = cbind(x_0,x_1,x_2,x_3,x_4)
  
  #psモデル
  #Data %*% t(beta)はXβ
  ps <- exp(Data %*% t(beta))/(1+exp(Data %*% t(beta)))
  
  #アウトカムモデル生成
  #delta = 2.5
  #zは介入,dataに生成したデータ全てをまとめる
  z <- rbinom(n,1,ps)
  y_out <- z*delta + Data %*% t(alpha) + matrix(rnorm(n,0,1),n,1)
  data <- cbind(Data,z,y_out)
  #--------------------------------------------------------------------#
  
  #傾向スコア推定
  res <- as.data.frame(data) %>%
    glm(z~x_0+x_1+x_2+x_3+x_4,data=.,family=binomial(link="logit"))
  ps_e <- res$fitted.values
  
  #Matching
  matching_out <- Match(Y= data[,7], Tr = data[,6], X =ps_e, 
                        caliper= 0.2,ties=FALSE,replace=FALSE)
  Mean[j,1] <- matching_out$est
  
  #Matching_Weight
  for (i in 1:n){
    w[i] = min(ps_e[i],1-ps_e[i])/(data[i,6]*ps_e[i] + (1-data[i,6])*(1-ps_e[i]))
  }
  Mean[j,2] = sum(w*data[,6]*data[,7])/sum(w*data[,6]) - sum(w*(1-data[,6])*data[,7])/sum(w*(1-data[,6]))
  
  boot_sd[j,] <- SD_boot(B = 20,d = data)
  
  
  cat("\r",j)
}
#----------------------------------------------------------------------------------------------------------#
#ボートストラップ標準誤差、信頼区間計算
CI = CI_normal(Mean = Mean, boot_sd=boot_sd)

#----------------------------------------------------------------------------------------------------------#


#結果

#Matching Weight
#bias
cat((mean(Mean[,2])-2.5)/2.5)
#MSE
cat(sum((Mean[,2]- 2.5)^2)/n)
#Emp SD
cat(sqrt(sum((Mean[,2]-mean(Mean[,2]))^2)/(n-1)))
#boot SD
cat(mean(boot_sd[,2]))
#CI
cat(CI[2])
#----------------------------#

#PS Matching
#bias
cat((mean(Mean[,1])-2.5)/2.5)
#MSE
cat(sum((Mean[,1]- 2.5)^2)/n)
#Emp SD
cat(sqrt(sum((Mean[,1]-mean(Mean[,1]))^2)/(n-1)))
#boot SD
cat(mean(boot_sd[,1]))
#CI
cat(CI[1])

