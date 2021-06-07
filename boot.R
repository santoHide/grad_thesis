  #ブートストラップ標準誤差
  SD_boot <- function(B=100,n=1000,d) {
    w_b <- rep(0,n)
    Mean_k_b <- rep(0,B)
    Mean_MW_b <- rep(0,B)
    b <-0
    if(is.null(d) ==FALSE) {
      for (b in 1:B){
        set.seed(b)
        n_b <- sample(nrow(d),n,replace=TRUE) 
        data_b <- as.data.frame(d[n_b,])
        #誤特定:X_2を取り除く
        res <-suppressWarnings(glm(formula=z~x_1+x_2+x_3+x_4,data= data_b,family=binomial(link="logit")))
        ps_e_b <- res$fitted.values
        #PSmatching
        m_out_b <-  Match(Y= data_b[,7], Tr = data_b[,6], X =ps_e_b, caliper= 0.2,ties=FALSE,replace=FALSE)
        Mean_k_b[b]<- m_out_b$est
        #Matching_Weight
        for (i in 1:n){
          w_b[i] <- min(ps_e_b[i],1-ps_e_b[i])/(data_b[i,6]*ps_e_b[i] + (1-data_b[i,6])*(1-ps_e_b[i]))
        }
        Mean_MW_b[b] <- sum(w_b*data_b[,6]*data_b[,7])/sum(w_b*data_b[,6]) - 
          sum(w_b*(1-data_b[,6])*data_b[,7])/sum(w_b*(1-data_b[,6]))
        
      }
    }
    r <- rep(0,4)
    r[1] <- mean(Mean_k_b)
    r[2] <- mean(Mean_MW_b)
    r[3] <- sqrt(sum((Mean_k_b - mean(Mean_k_b))^2)/(B-1))
    r[4] <- sqrt(sum((Mean_MW_b - mean(Mean_MW_b))^2)/(B-1))
    return(r)
  }




