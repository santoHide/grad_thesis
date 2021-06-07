CI_normal <- function(Mean= NULL, boot_sd = NULL) {
  CI <- rep(0,2)
  if(is.null(Mean) == FALSE && is.null(boot_sd) == FALSE){
    i <- 0
    for  (i in 1:dim(Mean)[1]){
      print(i)
      if(boot_sd[i,1] - 1.96*boot_sd[i,3] <= Mean[i,1] && Mean[i,1] <= boot_sd[i,1] + 1.96*boot_sd[i,3]){
        CI[1] = CI[1] + 1
      }
      if(boot_sd[i,2] - 1.96*boot_sd[i,4] <= Mean[i,2] && Mean[i,2] <= boot_sd[i,2] + 1.96*boot_sd[i,4]){
        CI[2] = CI[2] + 1
      }
    }
  }
  else{
     #エラー処理
    stop("ERROR\n")
  }
  return(CI)
}