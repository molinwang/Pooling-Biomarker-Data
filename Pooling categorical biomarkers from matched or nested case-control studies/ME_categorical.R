################################################################
####  The R function ME_categorical()                       ####
####  implements the methods proposed by Wu et al. for      ####
####  pooling categorical biomarker measurements            ####
####  from multi-center matched/nested case-control studies ####
################################################################

library(nnet)
#################################################################
################## Useful intermediate functions ################
#################################################################
lik_nk_fun_strata <- function(p, n, beta_mle0, cali_prob){
  #### n corrsponds to the n-th case in the study
  ####    (or the n-th stratum)
  #### p is the p_1, \ldots, p_{M_{j(s)}+1} indicator, with
  #### first element be the case and the remaining one be the ctrls
  cali_prob.temp <- cali_prob[cali_prob[,"strat"] == n,]
  num_match      <- nrow(cali_prob.temp) - 1 #### -1: the case
  p              <- unlist(p)
  #### beta_mle0 are the parameters to be estimated
  beta_mle       <- c(0,beta_mle0)
  #### likelihood in one stratum for one combination of p_1,\ldots, p_{M_{j(s)}+1} 
  (exp(beta_mle[p[1]+1])/sum(sapply(0:num_match,function(x)exp(beta_mle[p[x+1]+1]))))*
    prod(sapply(0:num_match,function(x)cali_prob.temp[x+1,p[x+1]+1]))
}     
lik_k_fun         <- function(beta_mle0, cali_prob){
  ### for each stratum, we first need summing up
  ### all possible combination of p_1,\ldots,p_{M_{j(s)+1}}
  #### P is the total number of categoriees
  log_lik_k <- NULL
  ncase     <- length(unique(cali_prob[,"strat"]))
  P         <- ncol(cali_prob) - 1 # calib_prob contains probs for each category + stratum number
  for(n in 1:ncase){
    num_match <- nrow(cali_prob[cali_prob[,"strat"]==n,]) - 1 ### -1 for the case
    p_eval <- expand.grid(replicate( num_match + 1, 0:(P-1), simplify=FALSE))
    log_lik_k[n] <- log(sum(sapply(split(p_eval, seq(nrow(p_eval))),lik_nk_fun_strata,
                                   n = n, beta_mle0 = beta_mle0,
                                   cali_prob = cali_prob)))
    
  }
  return(sum(log_lik_k))
}
lik_fun           <- function(beta_mle0, cali_prob_lab){
  log_lik <- NULL
  for (i in 1:length(cali_prob_lab)){
    log_lik[i]<-lik_k_fun(beta_mle0 = beta_mle0,
                          cali_prob = cali_prob_lab[[i]])
  }
  return(-sum(log_lik))
}  
psi_a_q_sj        <- function(reference, param_est, cali_prob){
  # Notice here the category starts from 0
  # reference: the q-th calibration subset
  # param_est: calibration model estimation
  # cali_prob is used for getting the number of strata
  ncase <- length(unique(cali_prob[,"strat"]))
  a <- matrix(0, nrow = nrow(param_est),
              ncol = ncase)
  for(i in 1:nrow(param_est)){
    temp <- ((reference$X_cate==i) - 
               exp(param_est[i,1] + param_est[i,2]*reference$W)/
               (1+ rowSums(sapply(1:nrow(param_est), 
                                  function(x)exp(param_est[x,1] + param_est[x,2]*reference$W)))))
    ### for 1:2 or more, more than one controls may be sampled
    temp2 <- sapply( reference$strat, function(x)sum(temp[reference$strat==x]))
    a[i,reference$strat] <- temp2
  }
  return(a)
}
psi_b_q_sj        <- function(reference, param_est, cali_prob){
  # reference: the q-th calibration subset
  # param_est: calibration model estimation
  ncase <- length(unique(cali_prob[,"strat"]))
  b <- matrix(0, nrow = nrow(param_est),
              ncol = ncase)
  for(i in 1:nrow(param_est)){
    temp <-(reference$W)*((reference$X_cate==i) - 
                            exp(param_est[i,1] + param_est[i,2]*reference$W)/
                            (1+ rowSums(sapply(1:nrow(param_est), 
                                               function(x)exp(param_est[x,1] + param_est[x,2]*reference$W)))))
    ### for 1:2 or more, more than one controls may be sampled
    temp2 <- sapply( reference$strat, function(x)sum(temp[reference$strat==x]))
    b[i,reference$strat] <- temp2
  }
  return(b)
}
psi_beta_x_sj_pi  <- function(p, n, beta_estimate, cali_prob){
  #### n corrsponds to the n-th case in the study
  ####    (or the n-th stratum)
  #### p is the p_2, \ldots, p_{M_{j(s)}+1} indicator, with
  #### first element be the case and the remaining one be the ctrls
  cali_prob.temp <- cali_prob[cali_prob[, "strat"] == n, ]
  num_match      <- nrow(cali_prob.temp) - 1
  #### notice category starts from 0
  p        <- unlist(p)
  p.case   <- p[1]  
  p.ctrl   <- p[-1] 
  #### beta_mle are the estimated parameters for beta_X 
  beta_mle <- matrix(beta_estimate, nrow = length(beta_estimate),
                     ncol = 1)
  #### B matrix for case
  B.case   <- matrix(0, nrow = nrow(beta_mle), ncol = 1)
  if(p.case!=0){
    B.case[p.case,] <- 1
  }
  B.ctrl   <- matrix(0, nrow = nrow(beta_mle), ncol = num_match)
  for(i in 1:num_match){
    if(p.ctrl[i]!=0){
      B.ctrl[p.ctrl[i],i] <- 1
    }
  }
  B.cbind  <- cbind(B.case, B.ctrl)
  
  psi_x    <- matrix(NA, nrow = length(beta_estimate), ncol = 1)
  useful_q <- sum(sapply(1:(num_match+1), function(x)exp(t(beta_mle)%*% B.cbind[,x])))
  #### likelihood in one stratum for one combination of p_1,\ldots, p_{M_{j(s)}+1} 
  
  for(p.prime in 1:length(beta_estimate)){
    numerator.1 <- c(B.case[p.prime,])*exp(t(beta_mle)%*% B.case)*useful_q
    
    numerator.2 <- c(exp(t(beta_mle)%*%B.case))*
      sum(sapply(1:(num_match+1), 
                 function(x)c(B.cbind[p.prime,x])*c(exp(t(beta_mle)%*%B.cbind[,x]))))
    
    numerator.3 <- prod(sapply(0:num_match, function(x)cali_prob.temp[x+1, p[x+1]+1]))
    
    denominator <- useful_q^2
    
    psi_x[p.prime,] <- ((numerator.1-numerator.2)*numerator.3)/denominator
  }
  return(psi_x)
}
psi_beta_x_sj     <- function(beta_estimate, cali_prob){
  ncase <- length(unique(cali_prob[,"strat"]))
  P         <- ncol(cali_prob) - 1 # calib_prob contains probs for each category + stratum number
  
  ### for each stratum, we need summing up
  ### all possible combination of p_1,\ldots,p_{M_{j(s)+1}}
  psi_beta_x_sj <- matrix(nrow = length(beta_estimate), ncol = ncase)
  for(n in 1:ncase){
    cali_prob.temp <- cali_prob[cali_prob[, "strat"] == n, ]
    num_match      <- nrow(cali_prob.temp) - 1
    p_eval <- expand.grid(replicate( num_match + 1, 0:(P-1), simplify=FALSE))
    L_sj   <- sum(sapply(split(p_eval, seq(nrow(p_eval))), lik_nk_fun_strata,
                         n = n, beta_mle0 = beta_estimate,
                         cali_prob = cali_prob))
    psi_beta_x_sj[,n]    <- matrix(rowSums(sapply(split(p_eval, seq(nrow(p_eval))),
                                                  psi_beta_x_sj_pi,
                                                  n = n, beta_estimate = beta_estimate,
                                                  cali_prob = cali_prob)), ncol = 1)/L_sj
  }
  return(psi_beta_x_sj)
}
psi_double_a_q_s  <- function(reference, param_est){
  # Notice here the category starts from 0
  # reference: the q-th calibration subset
  # param_est: calibration model estimation
  double_a <- matrix(NA, nrow = nrow(param_est), 
                     ncol = nrow(param_est))
  
  for(p.prime in 1:nrow(param_est)){
    for(p.star in 1:nrow(param_est)){
      if(p.prime == p.star){
        numerator.1.1 <- exp(param_est[p.prime,1] + param_est[p.prime,2]*reference$W)
        numerator.1.2 <- 1 + rowSums(sapply(1:nrow(param_est), 
                                            function(x)exp(param_est[x,1] + param_est[x,2]*reference$W)))
        numerator.2.sq   <- numerator.1.1
        #num           <- numerator.1.1*numerator.1.2 - numerator.2
        #denom         <- numerator.1.2^2
        frac1         <- numerator.1.1/numerator.1.2 
        frac2         <- (numerator.2.sq/numerator.1.2)^2
        double_a[p.prime, p.star] <- sum(-(frac1 - frac2))
      } else{
        numerator.1.1 <- exp(param_est[p.prime,1] + param_est[p.prime,2]*reference$W)
        numerator.1.2 <- exp(param_est[p.star,1]  + param_est[p.star,2]*reference$W)
        #num           <- numerator.1.1*numerator.1.2
        denom.sqrt    <- (1 + rowSums(sapply(1:nrow(param_est), 
                                             function(x)exp(param_est[x,1] + param_est[x,2]*reference$W))))
        frac1  <- numerator.1.1/denom.sqrt
        frac2  <- numerator.1.2/denom.sqrt
        
        double_a[p.prime, p.star] <- sum(frac1*frac2)
      }
    }
  }
  return(double_a)
}
psi_double_b_q_s  <- function(reference, param_est){
  # Notice here the category starts from 0
  # reference: the q-th calibration subset
  # param_est: calibration model estimation
  double_b <- matrix(NA, nrow = nrow(param_est), 
                     ncol = nrow(param_est))
  
  for(p.prime in 1:nrow(param_est)){
    for(p.star in 1:nrow(param_est)){
      if(p.prime == p.star){
        numerator.1.1 <- exp(param_est[p.prime,1] + param_est[p.prime,2]*reference$W)
        numerator.1.2 <- 1 + rowSums(sapply(1:nrow(param_est), 
                                            function(x)exp(param_est[x,1] + param_est[x,2]*reference$W)))
        numerator.2.sq   <- numerator.1.1
        #num           <- (reference$W^2)*(numerator.1.1*numerator.1.2 - numerator.2)
        # denom         <- numerator.1.2^2
        frac1         <- numerator.1.1/numerator.1.2 
        frac2         <- (numerator.2.sq/numerator.1.2)^2
        double_b[p.prime, p.star] <- sum(-(reference$W^2)*(frac1 - frac2))
      } else{
        numerator.1.1 <- exp(param_est[p.prime,1] + param_est[p.prime,2]*reference$W)
        numerator.1.2 <- exp(param_est[p.star,1]  + param_est[p.star,2]*reference$W)
        #num           <- (reference$W^2)*numerator.1.1*numerator.1.2
        denom.sqrt         <- (1 + rowSums(sapply(1:nrow(param_est), 
                                                  function(x)exp(param_est[x,1] + param_est[x,2]*reference$W))))
        frac1  <- numerator.1.1/denom.sqrt
        frac2  <- numerator.1.2/denom.sqrt
        double_b[p.prime, p.star] <- sum((reference$W^2)*frac1*frac2)
      }
    }
  }
  return(double_b)
}
psi_a_b_q_s       <- function(reference, param_est){
  # Notice here the category starts from 0
  # reference: the q-th calibration subset
  # param_est: calibration model estimation
  a_b <- matrix(NA, nrow = nrow(param_est), 
                ncol = nrow(param_est))
  ### In the simulation, we only have one b for each category
  for(p.prime in 1:nrow(param_est)){
    for(p.star in 1:nrow(param_est)){
      if(p.prime == p.star){
        numerator.1.1 <- reference$W*exp(param_est[p.prime,1] + param_est[p.prime,2]*reference$W)
        numerator.1.2 <- 1 + rowSums(sapply(1:nrow(param_est), 
                                            function(x)exp(param_est[x,1] + param_est[x,2]*reference$W)))
        #numerator.2   <- reference$W*(numerator.1.1^2)
        #num           <- numerator.1.1*numerator.1.2 - numerator.2
        #denom         <- numerator.1.2^2
        #num                <- numerator.1.1*numerator.1.2 - numerator.2
        denom.sqrt         <- numerator.1.2
        frac1              <- numerator.1.1/denom.sqrt
        frac2              <- reference$W*(numerator.1.1/denom.sqrt)^2
        a_b[p.prime, p.star] <- sum(-(frac1-frac2))
      } else{
        numerator.1.1 <- exp(param_est[p.prime,1] + param_est[p.prime,2]*reference$W)
        numerator.1.2 <- exp(param_est[p.star,1]  + param_est[p.star,2]*reference$W)
        #num           <- reference$W*numerator.1.1*numerator.1.2
        denom.sqrt         <- (1 + rowSums(sapply(1:nrow(param_est), 
                                                  function(x)exp(param_est[x,1] + param_est[x,2]*reference$W))))
        frac1         <- numerator.1.1/denom.sqrt
        frac2         <- numerator.1.2/denom.sqrt
        a_b[p.prime, p.star] <- sum(reference$W*frac1*frac2)
      }
    }
  }
  return(a_b)
}
mini.func.1_full  <- function(p, n, beta_estimate, cali_prob, param_est,
                              data_case, data_ctrl,
                              index.a){
  #### n corrsponds to the n-th case in the study
  ####    (or the n-th stratum)
  #### p is the p_2, \ldots, p_{M_{j(s)}+1} indicator, with
  #### first element be the case and the remaining one be the ctrls
  cali_prob.temp <- cali_prob[cali_prob[, "strat"] == n, ]
  num_match      <- nrow(cali_prob.temp) - 1
  #### notice category starts from 0
  p        <- unlist(p)
  p.case   <- p[1]  
  p.ctrl   <- p[-1] 
  #### beta_mle are the estimated parameters for beta_X 
  beta_mle <- matrix(beta_estimate, nrow = length(beta_estimate),
                     ncol = 1)
  #### B matrix for case
  B.case   <- matrix(0, nrow = nrow(beta_mle), ncol = 1)
  if(p.case!=0){
    B.case[p.case,] <- 1
  }
  B.ctrl   <- matrix(0, nrow = nrow(beta_mle), ncol = num_match)
  for(i in 1:num_match){
    if(p.ctrl[i]!=0){
      B.ctrl[p.ctrl[i],i] <- 1
    }
  }
  B.cbind  <- cbind(B.case, B.ctrl)
  
  temp_case <- data_case[data_case[,"strat"] == n,]
  temp_ctrl <- data_ctrl[data_ctrl[,"strat"] == n,]
  temp_cali_prob <- cali_prob.temp
  W.case    <- temp_case$W
  W.ctrl    <- temp_ctrl$W
  W.combined <- c(W.case, W.ctrl)
  
  numerator.3 <- 0
  for(i in 1:(num_match+1)){
    sub_denom <- 1 + sum(sapply(1:nrow(param_est), 
                                function(x)exp(param_est[x,1]+param_est[x,2]*W.combined[i])))
    sub_num   <- exp(param_est[index.a,1] + param_est[index.a,2]*W.combined[i])
    
    ################### p starts from 0 not 1!
    if(p[i]==0){
      prt.1  <- (-sub_num/sub_denom)/(sub_denom)
    } else if(p[i]==index.a){
      #prt.1  <- ( sub_num*sub_denom - sub_num^2 )/sub_denom^2
      prt.1  <- ( sub_num/sub_denom - (sub_num/sub_denom)^2 )
    } else if( (p[i]>0) & (p[i]!=index.a) ) {
      sub_num_2 <- exp(param_est[p[i],1] + param_est[p[i],2]*W.combined[i])
      # prt.1  <- (-sub_num*sub_num_2)/sub_denom^2
      prt.1  <- -(sub_num/sub_denom)*(sub_num_2/sub_denom)
    }
    
    prt.2  <-  prod(sapply(c(1:(num_match+1))[-i],
                           function(x)temp_cali_prob[x, p[x]+1])) 
    
    numerator.3 <- numerator.3 + prt.1*prt.2
  }
  prt.3  <- exp(t(beta_mle)%*%B.case)/
    (sum(sapply(1:(num_match+1), function(x)exp(t(beta_mle)%*% B.cbind[,x]))))
  numerator.4 <- numerator.3 * prt.3
  return(numerator.4)
}
mini.func.1_internal        <- function(p, n, beta_estimate, cali_prob, param_est,
                                        data_case, data_ctrl,
                                        index.a){
  #### n corrsponds to the n-th case in the study
  ####    (or the n-th stratum)
  #### p is the p_2, \ldots, p_{M_{j(s)}+1} indicator, with
  #### first element be the case and the remaining one be the ctrls
  cali_prob.temp <- cali_prob[cali_prob[, "strat"] == n, ]
  num_match      <- nrow(cali_prob.temp) - 1
  #### notice category starts from 0
  p        <- unlist(p)
  p.case   <- p[1]  
  p.ctrl   <- p[-1] 
  #### beta_mle are the estimated parameters for beta_X 
  beta_mle <- matrix(beta_estimate, nrow = length(beta_estimate),
                     ncol = 1)
  #### B matrix for case
  B.case   <- matrix(0, nrow = nrow(beta_mle), ncol = 1)
  if(p.case!=0){
    B.case[p.case,] <- 1
  }
  B.ctrl   <- matrix(0, nrow = nrow(beta_mle), ncol = num_match)
  for(i in 1:num_match){
    if(p.ctrl[i]!=0){
      B.ctrl[p.ctrl[i],i] <- 1
    }
  }
  B.cbind  <- cbind(B.case, B.ctrl)
  
  temp_case <- data_case[data_case[,"strat"] == n,]
  temp_ctrl <- data_ctrl[data_ctrl[,"strat"] == n,]
  temp_cali_prob <- cali_prob.temp
  W.case    <- temp_case$W
  W.ctrl    <- temp_ctrl$W
  W.combined <- c(W.case, W.ctrl)
  C.combined <- c(temp_case$C, temp_ctrl$C)
  
  numerator.3 <- 0
  for(i in 1:(num_match+1)){
    sub_denom <- 1 + sum(sapply(1:nrow(param_est), 
                                function(x)exp(param_est[x,1]+param_est[x,2]*W.combined[i])))
    sub_num   <- exp(param_est[index.a,1] + param_est[index.a,2]*W.combined[i])
    
    ################### p starts from 0 not 1!
    if(C.combined[i] == FALSE){
      if(p[i]==0){
        prt.1  <- (-sub_num/sub_denom)/(sub_denom)
      } else if(p[i]==index.a){
        prt.1  <- (sub_num/sub_denom) - (sub_num/sub_denom)^2 
      } else if( (p[i]>0) & (p[i]!=index.a) ) {
        sub_num_2 <- exp(param_est[p[i],1] + param_est[p[i],2]*W.combined[i])
        prt.1  <- -(sub_num/sub_denom)*(sub_num_2/sub_denom)
      }
    } else{
      prt.1 <- 0
    }
    
    prt.2  <-  prod(sapply(c(1:(num_match+1))[-i],
                           function(x)temp_cali_prob[x, p[x]+1])) 
    
    numerator.3 <- numerator.3 + prt.1*prt.2
  }
  prt.3  <- exp(t(beta_mle)%*%B.case)/
    (sum(sapply(1:(num_match+1), function(x)exp(t(beta_mle)%*% B.cbind[,x]))))
  numerator.4 <- numerator.3 * prt.3
  return(numerator.4)
}
psi_beta_x_a_sj_pi_full     <- function(p, n, beta_estimate, cali_prob, param_est,
                                        data_case, data_ctrl){
  cali_prob.temp <- cali_prob[cali_prob[, "strat"] == n, ]
  num_match      <- nrow(cali_prob.temp) - 1
  P              <- ncol(cali_prob) - 1 # calib_prob contains probs for each category + stratum number
  
  #### n corrsponds to the n-th case in the study
  ####    (or the n-th stratum)
  #### p is the p_2, \ldots, p_{M_{j(s)}+1} indicator, with
  #### first element be the case and the remaining one be the ctrls
  p_eval <- expand.grid(replicate( num_match + 1, 0:(P-1), simplify=FALSE))
  L_sj   <- sum(sapply(split(p_eval, seq(nrow(p_eval))), lik_nk_fun_strata,
                       n = n, beta_mle0 = beta_estimate,
                       cali_prob = cali_prob))
  
  #### notice category starts from 0
  p        <- unlist(p)
  p.case   <- p[1]  
  p.ctrl   <- p[-1] 
  #### beta_mle are the estimated parameters for beta_X 
  beta_mle <- matrix(beta_estimate, nrow = length(beta_estimate),
                     ncol = 1)
  #### B matrix for case
  B.case   <- matrix(0, nrow = nrow(beta_mle), ncol = 1)
  if(p.case!=0){
    B.case[p.case,] <- 1
  }
  B.ctrl   <- matrix(0, nrow = nrow(beta_mle), ncol = num_match)
  for(i in 1:num_match){
    if(p.ctrl[i]!=0){
      B.ctrl[p.ctrl[i],i] <- 1
    }
  }
  B.cbind  <- cbind(B.case, B.ctrl)
  
  temp_case <- data_case[data_case[,"strat"] == n,]
  temp_ctrl <- data_ctrl[data_ctrl[,"strat"] == n,]
  temp_cali_prob <- cali_prob.temp
  W.case    <- temp_case$W
  W.ctrl    <- temp_ctrl$W
  W.combined <- c(W.case, W.ctrl)
  #### object to store the second order derivative
  x_a <- matrix(NA, nrow = length(beta_estimate),
                ncol = nrow(param_est))
  for(index.a in 1:nrow(param_est)){
    numerator.1 <- 0
    #### Calculate sum_{i=1}^{M_{j(s)+1}}
    for(i in 1:(num_match+1)){
      sub_denom <- 1 + sum(sapply(1:nrow(param_est), 
                                  function(x)exp(param_est[x,1]+param_est[x,2]*W.combined[i])))
      sub_num   <- exp(param_est[index.a,1] + param_est[index.a,2]*W.combined[i])
      
      ################### p starts from 0 not 1!
      if(p[i]==0){
        prt.1  <- (-sub_num/sub_denom)/(sub_denom)
      } else if(p[i]==index.a){
        prt.1  <- ( sub_num/sub_denom) - (sub_num/sub_denom)^2 
      } else if( (p[i]>0) & (p[i]!=index.a) ) {
        sub_num_2 <- exp(param_est[p[i],1] + param_est[p[i],2]*W.combined[i])
        prt.1     <- -(sub_num/sub_denom)*(sub_num_2/sub_denom)
      }
      
      prt.2  <- prod(sapply(c(1:(num_match+1))[-i],
                            function(x)temp_cali_prob[x, p[x]+1])) 
      
      numerator.1 <- numerator.1 + prt.1*prt.2
    }
    numerator.1    <- numerator.1 * L_sj
    
    numerator.2    <- prod(sapply(1:(num_match+1), 
                                  function(x)temp_cali_prob[x, p[x]+1])) 
    
    numerator.3 <- sum(sapply( split(p_eval, seq(nrow(p_eval))), mini.func.1_full,
                               n=n, beta_estimate=beta_estimate, 
                               cali_prob = cali_prob, param_est = param_est,
                               data_case = data_case, data_ctrl = data_ctrl,
                               #ncase = ncase, num_match = num_match,
                               index.a = index.a))
    
    Numerator_1   <- numerator.1-numerator.2*numerator.3
    
    Numerator_2.1 <- B.case%*%(exp( t(beta_mle)%*% B.case)*
                                 sum(sapply(1:(num_match+1), 
                                            function(x)exp(t(beta_mle)%*% B.cbind[,x]))))
    
    Numerator_2.2 <- c(exp(t(beta_mle)%*%B.case))*
      matrix(rowSums(sapply(1:(num_match+1), 
                            function(x)B.cbind[,x]*c(exp(t(beta_mle)%*%B.cbind[,x])))), ncol = 1)
    
    Denominator_2 <-  (sum(sapply(1:(num_match+1), 
                                  function(x)exp(t(beta_mle)%*% B.cbind[,x]))))^2
    
    Numerator_2   <- ((Numerator_2.1/Denominator_2)- (Numerator_2.2/Denominator_2))/Denominator_2
    
    Numerator_all <- c(Numerator_1)*Numerator_2
    
    final.value <- (Numerator_all/L_sj)/(L_sj)
    
    x_a[, index.a] <- final.value
  }
  
  return( x_a )
}
psi_beta_x_a_sj_pi_internal <- function(p, n, beta_estimate, cali_prob, param_est,
                                        data_case, data_ctrl){
  cali_prob.temp <- cali_prob[cali_prob[, "strat"] == n, ]
  num_match      <- nrow(cali_prob.temp) - 1
  P              <- ncol(cali_prob) - 1 # calib_prob contains probs for each category + stratum number
  
  #### n corrsponds to the n-th case in the study
  ####    (or the n-th stratum)
  #### p is the p_2, \ldots, p_{M_{j(s)}+1} indicator, with
  #### first element be the case and the remaining one be the ctrls
  p_eval <- expand.grid(replicate( num_match + 1, 0:(P-1), simplify=FALSE))
  L_sj   <- sum(sapply(split(p_eval, seq(nrow(p_eval))), lik_nk_fun_strata,
                       n = n, beta_mle0 = beta_estimate,
                       cali_prob = cali_prob))
  
  #### notice category starts from 0
  p        <- unlist(p)
  p.case   <- p[1]  
  p.ctrl   <- p[-1] 
  #### beta_mle are the estimated parameters for beta_X 
  beta_mle <- matrix(beta_estimate, nrow = length(beta_estimate),
                     ncol = 1)
  #### B matrix for case
  B.case   <- matrix(0, nrow = nrow(beta_mle), ncol = 1)
  if(p.case!=0){
    B.case[p.case,] <- 1
  }
  B.ctrl   <- matrix(0, nrow = nrow(beta_mle), ncol = num_match)
  for(i in 1:num_match){
    if(p.ctrl[i]!=0){
      B.ctrl[p.ctrl[i],i] <- 1
    }
  }
  B.cbind  <- cbind(B.case, B.ctrl)
  
  temp_case <- data_case[data_case[,"strat"] == n,]
  temp_ctrl <- data_ctrl[data_ctrl[,"strat"] == n,]
  temp_cali_prob <- cali_prob.temp
  W.case    <- temp_case$W
  W.ctrl    <- temp_ctrl$W
  W.combined <- c(W.case, W.ctrl)
  C.combined <- c(temp_case$C, temp_ctrl$C)
  
  #### object to store the second order derivative
  x_a <- matrix(NA, nrow = length(beta_estimate),
                ncol = nrow(param_est))
  for(index.a in 1:nrow(param_est)){
    numerator.1 <- 0
    #### Calculate sum_{i=1}^{M_{j(s)+1}}
    for(i in 1:(num_match+1)){
      sub_denom <- 1 + sum(sapply(1:nrow(param_est), 
                                  function(x)exp(param_est[x,1]+param_est[x,2]*W.combined[i])))
      sub_num   <- exp(param_est[index.a,1] + param_est[index.a,2]*W.combined[i])
      
      ################### p starts from 0 not 1!
      if(C.combined[i]==FALSE){
        if(p[i]==0){
          prt.1  <- (-sub_num/sub_denom)/sub_denom
        } else if(p[i]==index.a){
          prt.1  <- ( sub_num/sub_denom) - (sub_num/sub_denom)^2 
        } else if( (p[i]>0) & (p[i]!=index.a) ) {
          sub_num_2 <- exp(param_est[p[i],1] + param_est[p[i],2]*W.combined[i])
          prt.1     <- (-sub_num/sub_denom)*(sub_num_2/sub_denom)
        }
      } else{
        prt.1 <- 0
      }
      prt.2  <- prod(sapply(c(1:(num_match+1))[-i],
                            function(x)temp_cali_prob[x, p[x]+1])) 
      
      numerator.1 <- numerator.1 + prt.1*prt.2
    }
    numerator.1    <- numerator.1 * L_sj
    
    numerator.2    <- prod(sapply(1:(num_match+1), 
                                  function(x)temp_cali_prob[x, p[x]+1])) 
    
    numerator.3 <- sum(sapply( split(p_eval, seq(nrow(p_eval))), mini.func.1_internal,
                               n=n, beta_estimate=beta_estimate, 
                               cali_prob = cali_prob, param_est = param_est,
                               data_case = data_case, data_ctrl = data_ctrl,
                               #ncase = ncase, num_match = num_match,
                               index.a = index.a))
    
    Numerator_1   <- numerator.1-numerator.2*numerator.3
    
    Numerator_2.1 <- B.case%*%(exp( t(beta_mle)%*% B.case)*
                                 sum(sapply(1:(num_match+1), 
                                            function(x)exp(t(beta_mle)%*% B.cbind[,x]))))
    
    Numerator_2.2 <- c(exp(t(beta_mle)%*%B.case))*
      matrix(rowSums(sapply(1:(num_match+1), 
                            function(x)B.cbind[,x]*c(exp(t(beta_mle)%*%B.cbind[,x])))), ncol = 1)
    
    Denominator_2 <-  (sum(sapply(1:(num_match+1), 
                                  function(x)exp(t(beta_mle)%*% B.cbind[,x]))))^2
    
    Numerator_2   <- ((Numerator_2.1/Denominator_2) - (Numerator_2.2/Denominator_2))/Denominator_2
    
    Numerator_all <- c(Numerator_1)*Numerator_2
    
    final.value <- (Numerator_all/L_sj)/(L_sj)
    
    x_a[, index.a] <- final.value
  }
  
  return( x_a )
}
psi_beta_x_a_sj             <- function(beta_estimate, cali_prob, param_est,
                                        data_case, data_ctrl, 
                                        calib_method){
  ncase       <- length(unique(cali_prob[,"strat"]))
  beta_x_a_sj <- list()
  P           <- ncol(cali_prob) - 1 # calib_prob contains probs for each category + stratum number
  
  for(n in 1:ncase){
    cali_prob.temp <- cali_prob[cali_prob[, "strat"] == n, ]
    num_match      <- nrow(cali_prob.temp) - 1
    p_eval <- expand.grid(replicate( num_match + 1, 0:(P-1), simplify=FALSE))
    #### sapply returns a vector by column instead of a matrix 
    if(calib_method=="Full"){
      beta_x_a_sj.temp <- rowSums(sapply(split(p_eval, seq(nrow(p_eval))), psi_beta_x_a_sj_pi_full,
                                         n=n, beta_estimate=beta_estimate, 
                                         cali_prob = cali_prob, param_est = param_est,
                                         data_case = data_case, data_ctrl = data_ctrl))
    } else if(calib_method=="Internalized"){
      beta_x_a_sj.temp <- rowSums(sapply(split(p_eval, seq(nrow(p_eval))), psi_beta_x_a_sj_pi_internal,
                                         n=n, beta_estimate=beta_estimate, 
                                         cali_prob = cali_prob, param_est = param_est,
                                         data_case = data_case, data_ctrl = data_ctrl))
    }
    beta_x_a_sj[[n]] <- matrix(beta_x_a_sj.temp, nrow = length(beta_estimate), 
                               ncol = nrow(param_est), byrow = FALSE)
  }
  return(beta_x_a_sj)
}
mini.func.2_full            <- function(p, n, beta_estimate, cali_prob, param_est,
                                        data_case, data_ctrl, 
                                        index.b){
  #### n corrsponds to the n-th case in the study
  ####    (or the n-th stratum)
  #### p is the p_2, \ldots, p_{M_{j(s)}+1} indicator, with
  #### first element be the case and the remaining one be the ctrls
  cali_prob.temp <- cali_prob[cali_prob[, "strat"] == n, ]
  num_match      <- nrow(cali_prob.temp) - 1
  #### notice category starts from 0
  p        <- unlist(p)
  p.case   <- p[1]  
  p.ctrl   <- p[-1] 
  #### beta_mle are the estimated parameters for beta_X 
  beta_mle <- matrix(beta_estimate, nrow = length(beta_estimate),
                     ncol = 1)
  #### B matrix for case
  B.case   <- matrix(0, nrow = nrow(beta_mle), ncol = 1)
  if(p.case!=0){
    B.case[p.case,] <- 1
  }
  B.ctrl   <- matrix(0, nrow = nrow(beta_mle), ncol = num_match)
  for(i in 1:num_match){
    if(p.ctrl[i]!=0){
      B.ctrl[p.ctrl[i],i] <- 1
    }
  }
  B.cbind  <- cbind(B.case, B.ctrl)
  
  temp_case <- data_case[data_case[,"strat"] == n,]
  temp_ctrl <- data_ctrl[data_ctrl[,"strat"] == n,]
  temp_cali_prob <- cali_prob.temp
  W.case    <- temp_case$W
  W.ctrl    <- temp_ctrl$W
  W.combined <- c(W.case, W.ctrl)
  
  numerator.3 <- 0
  for(i in 1:(num_match+1)){
    sub_denom <- 1 + sum(sapply(1:nrow(param_est), 
                                function(x)exp(param_est[x,1]+param_est[x,2]*W.combined[i])))
    sub_num   <- exp(param_est[index.b,1] + param_est[index.b,2]*W.combined[i])
    
    ################### p starts from 0 not 1!
    if(p[i]==0){
      prt.1  <- (-W.combined[i]*(sub_num/sub_denom))/sub_denom
    } else if(p[i]==index.b){
      prt.1  <- ( W.combined[i]*(sub_num/sub_denom) - 
                    W.combined[i]*(sub_num/sub_denom)^2 )
    } else if( (p[i]>0) & (p[i]!=index.b) ) {
      sub_num_2 <- exp(param_est[p[i],1] + param_est[p[i],2]*W.combined[i])
      prt.1  <- (-W.combined[i]*(sub_num/sub_denom)*(sub_num_2/sub_denom))
    }
    
    prt.2  <-  prod(sapply(c(1:(num_match+1))[-i],
                           function(x)temp_cali_prob[x, p[x]+1])) 
    
    numerator.3 <- numerator.3 + prt.1*prt.2
  }
  prt.3  <- exp(t(beta_mle)%*%B.case)/
    (sum(sapply(1:(num_match+1), function(x)exp(t(beta_mle)%*% B.cbind[,x]))))
  numerator.3 <- numerator.3 * prt.3
  return(numerator.3)
}
mini.func.2_internal        <- function(p, n, beta_estimate, cali_prob, param_est,
                                        data_case, data_ctrl, 
                                        index.b){
  #### n corrsponds to the n-th case in the study
  ####    (or the n-th stratum)
  #### p is the p_2, \ldots, p_{M_{j(s)}+1} indicator, with
  #### first element be the case and the remaining one be the ctrls
  cali_prob.temp <- cali_prob[cali_prob[, "strat"] == n, ]
  num_match      <- nrow(cali_prob.temp) - 1
  #### notice category starts from 0
  p        <- unlist(p)
  p.case   <- p[1]  
  p.ctrl   <- p[-1] 
  #### beta_mle are the estimated parameters for beta_X 
  beta_mle <- matrix(beta_estimate, nrow = length(beta_estimate),
                     ncol = 1)
  #### B matrix for case
  B.case   <- matrix(0, nrow = nrow(beta_mle), ncol = 1)
  if(p.case!=0){
    B.case[p.case,] <- 1
  }
  B.ctrl   <- matrix(0, nrow = nrow(beta_mle), ncol = num_match)
  for(i in 1:num_match){
    if(p.ctrl[i]!=0){
      B.ctrl[p.ctrl[i],i] <- 1
    }
  }
  B.cbind  <- cbind(B.case, B.ctrl)
  
  temp_case <- data_case[data_case[,"strat"] == n,]
  temp_ctrl <- data_ctrl[data_ctrl[,"strat"] == n,]
  temp_cali_prob <- cali_prob.temp
  W.case    <- temp_case$W
  W.ctrl    <- temp_ctrl$W
  W.combined <- c(W.case, W.ctrl)
  C.combined <- c(temp_case$C, temp_ctrl$C)
  
  numerator.3 <- 0
  for(i in 1:(num_match+1)){
    sub_denom <- 1 + sum(sapply(1:nrow(param_est), 
                                function(x)exp(param_est[x,1]+param_est[x,2]*W.combined[i])))
    sub_num   <- exp(param_est[index.b,1] + param_est[index.b,2]*W.combined[i])
    
    ################### p starts from 0 not 1!
    if(C.combined[i]==FALSE){
      if(p[i]==0){
        prt.1  <- (-W.combined[i]*sub_num/sub_denom)/sub_denom
      } else if(p[i]==index.b){
        prt.1  <- ( W.combined[i]*(sub_num/sub_denom) - 
                      W.combined[i]*(sub_num/sub_denom)^2 )
      } else if( (p[i]>0) & (p[i]!=index.b) ) {
        sub_num_2 <- exp(param_est[p[i],1] + param_est[p[i],2]*W.combined[i])
        prt.1  <- (-W.combined[i]*(sub_num/sub_denom)*(sub_num_2/sub_denom))
      }
    } else{
      prt.1 <- 0
    }
    
    prt.2  <-  prod(sapply(c(1:(num_match+1))[-i],
                           function(x)temp_cali_prob[x, p[x]+1])) 
    
    numerator.3 <- numerator.3 + prt.1*prt.2
  }
  prt.3  <- exp(t(beta_mle)%*%B.case)/
    (sum(sapply(1:(num_match+1), function(x)exp(t(beta_mle)%*% B.cbind[,x]))))
  numerator.3 <- numerator.3 * prt.3
  return(numerator.3)
}
psi_beta_x_b_sj_pi_full     <- function(p, n, beta_estimate, cali_prob, param_est,
                                        data_case, data_ctrl){
  #### n corrsponds to the n-th case in the study
  ####    (or the n-th stratum)
  #### p is the p_2, \ldots, p_{M_{j(s)}+1} indicator, with
  #### first element be the case and the remaining one be the ctrls
  cali_prob.temp <- cali_prob[cali_prob[, "strat"] == n, ]
  num_match      <- nrow(cali_prob.temp) - 1
  P              <- ncol(cali_prob) - 1 # calib_prob contains probs for each category + stratum number
  
  p_eval <- expand.grid(replicate( num_match + 1, 0:(P-1), simplify=FALSE))
  L_sj   <- sum(sapply(split(p_eval, seq(nrow(p_eval))), lik_nk_fun_strata,
                       n = n, beta_mle0 = beta_estimate,
                       cali_prob = cali_prob))
  
  #### notice category starts from 0
  p        <- unlist(p)
  p.case   <- p[1]  
  p.ctrl   <- p[-1] 
  #### beta_mle are the estimated parameters for beta_X 
  beta_mle <- matrix(beta_estimate, nrow = length(beta_estimate),
                     ncol = 1)
  #### B matrix for case
  B.case   <- matrix(0, nrow = nrow(beta_mle), ncol = 1)
  if(p.case!=0){
    B.case[p.case,] <- 1
  }
  B.ctrl   <- matrix(0, nrow = nrow(beta_mle), ncol = num_match)
  for(i in 1:num_match){
    if(p.ctrl[i]!=0){
      B.ctrl[p.ctrl[i],i] <- 1
    }
  }
  B.cbind  <- cbind(B.case, B.ctrl)
  
  temp_case <- data_case[data_case[,"strat"] == n,]
  temp_ctrl <- data_ctrl[data_ctrl[,"strat"] == n,]
  temp_cali_prob <- cali_prob.temp
  W.case    <- temp_case$W
  W.ctrl    <- temp_ctrl$W
  W.combined <- c(W.case, W.ctrl)
  #### object to store the second order derivative
  x_b <- matrix(NA, nrow = length(beta_estimate),
                ncol = nrow(param_est))
  for(index.b in 1:nrow(param_est)){
    numerator.1 <- 0
    #### Calculate sum_{i=1}^{M_{j(s)+1}}
    for(i in 1:(num_match+1)){
      sub_denom <- 1 + sum(sapply(1:nrow(param_est), 
                                  function(x)exp(param_est[x,1]+param_est[x,2]*W.combined[i])))
      sub_num   <- exp(param_est[index.b,1] + param_est[index.b,2]*W.combined[i])
      
      ################### p starts from 0 not 1!
      if(p[i]==0){
        prt.1  <- (-W.combined[i]*sub_num/sub_denom)/(sub_denom)
      } else if(p[i]==index.b){
        prt.1  <- ( W.combined[i]*(sub_num/sub_denom) - 
                      W.combined[i]*(sub_num/sub_denom)^2 )
      } else if( (p[i]>0) & (p[i]!=index.b) ) {
        sub_num_2 <- exp(param_est[p[i],1] + param_est[p[i],2]*W.combined[i])
        prt.1  <- (-W.combined[i]*(sub_num/sub_denom)*(sub_num_2/sub_denom))
      }
      
      prt.2  <- prod(sapply(c(1:(num_match+1))[-i],
                            function(x)temp_cali_prob[x, p[x]+1])) 
      
      numerator.1 <- numerator.1 + prt.1*prt.2
    }
    numerator.1    <- numerator.1 * L_sj
    
    numerator.2    <- prod(sapply(1:(num_match+1), 
                                  function(x)temp_cali_prob[x, p[x]+1])) 
    
    numerator.3 <- sum(sapply( split(p_eval, seq(nrow(p_eval))), mini.func.2_full,
                               n=n, beta_estimate=beta_estimate, 
                               cali_prob = cali_prob, param_est = param_est,
                               data_case = data_case, data_ctrl = data_ctrl,
                               index.b = index.b))
    
    Numerator_1   <- numerator.1-numerator.2*numerator.3
    
    Numerator_2.1 <- B.case%*%(exp(t(beta_mle)%*%B.case)*
                                 sum(sapply(1:(num_match+1), function(x)exp(t(beta_mle)%*%B.cbind[,x]))))
    
    Numerator_2.2 <- c(exp(t(beta_mle)%*% B.case))*
      matrix(rowSums(sapply(1:(num_match+1), 
                            function(x)B.cbind[,x]*c(exp(t(beta_mle)%*% B.cbind[,x])))), ncol = 1)
    
    Denominator_2 <-  (sum(sapply(1:(num_match+1), 
                                  function(x)exp(t(beta_mle)%*% B.cbind[,x]))))^2
    
    Numerator_2   <- (Numerator_2.1/Denominator_2- Numerator_2.2/Denominator_2)/Denominator_2
    
    Numerator_all <- c(Numerator_1)*Numerator_2
    
    final.value <- (Numerator_all/L_sj)/L_sj
    
    x_b[, index.b] <- final.value
  }
  
  return( x_b )
}

psi_beta_x_b_sj_pi_internal <- function(p, n, beta_estimate, cali_prob, param_est,
                                        data_case, data_ctrl){
  #### n corrsponds to the n-th case in the study
  ####    (or the n-th stratum)
  #### p is the p_2, \ldots, p_{M_{j(s)}+1} indicator, with
  #### first element be the case and the remaining one be the ctrls
  cali_prob.temp <- cali_prob[cali_prob[, "strat"] == n, ]
  num_match      <- nrow(cali_prob.temp) - 1
  P              <- ncol(cali_prob) - 1 # calib_prob contains probs for each category + stratum number
  
  p_eval <- expand.grid(replicate( num_match + 1, 0:(P-1), simplify=FALSE))
  L_sj   <- sum(sapply(split(p_eval, seq(nrow(p_eval))), lik_nk_fun_strata,
                       n = n, beta_mle0 = beta_estimate,
                       cali_prob = cali_prob))
  
  #### notice category starts from 0
  p        <- unlist(p)
  p.case   <- p[1]  
  p.ctrl   <- p[-1] 
  #### beta_mle are the estimated parameters for beta_X 
  beta_mle <- matrix(beta_estimate, nrow = length(beta_estimate),
                     ncol = 1)
  #### B matrix for case
  B.case   <- matrix(0, nrow = nrow(beta_mle), ncol = 1)
  if(p.case!=0){
    B.case[p.case,] <- 1
  }
  B.ctrl   <- matrix(0, nrow = nrow(beta_mle), ncol = num_match)
  for(i in 1:num_match){
    if(p.ctrl[i]!=0){
      B.ctrl[p.ctrl[i],i] <- 1
    }
  }
  B.cbind  <- cbind(B.case, B.ctrl)
  
  temp_case <- data_case[data_case[,"strat"] == n,]
  temp_ctrl <- data_ctrl[data_ctrl[,"strat"] == n,]
  temp_cali_prob <- cali_prob.temp
  W.case    <- temp_case$W
  W.ctrl    <- temp_ctrl$W
  W.combined <- c(W.case, W.ctrl)
  C.combined <- c(temp_case$C, temp_ctrl$C)
  #### object to store the second order derivative
  x_b <- matrix(NA, nrow = length(beta_estimate),
                ncol = nrow(param_est))
  for(index.b in 1:nrow(param_est)){
    numerator.1 <- 0
    #### Calculate sum_{i=1}^{M_{j(s)+1}}
    for(i in 1:(num_match+1)){
      sub_denom <- 1 + sum(sapply(1:nrow(param_est), 
                                  function(x)exp(param_est[x,1]+param_est[x,2]*W.combined[i])))
      sub_num   <- exp(param_est[index.b,1] + param_est[index.b,2]*W.combined[i])
      
      ################### p starts from 0 not 1!
      if(C.combined[i]==FALSE){
        if(p[i]==0){
          prt.1  <- (-W.combined[i]*sub_num/sub_denom)/sub_denom
        } else if(p[i]==index.b){
          prt.1  <- ( W.combined[i]*(sub_num/sub_denom) - 
                        W.combined[i]*(sub_num/sub_denom)^2 )
        } else if( (p[i]>0) & (p[i]!=index.b) ) {
          sub_num_2 <- exp(param_est[p[i],1] + param_est[p[i],2]*W.combined[i])
          prt.1  <- (-W.combined[i]*(sub_num/sub_denom)*(sub_num_2/sub_denom))
        }
      } else{
        prt.1 <- 0
      }
      
      prt.2  <- prod(sapply(c(1:(num_match+1))[-i],
                            function(x)temp_cali_prob[x, p[x]+1])) 
      
      numerator.1 <- numerator.1 + prt.1*prt.2
    }
    numerator.1    <- numerator.1 * L_sj
    
    numerator.2    <- prod(sapply(1:(num_match+1), 
                                  function(x)temp_cali_prob[x, p[x]+1])) 
    
    numerator.3 <- sum(sapply( split(p_eval, seq(nrow(p_eval))), mini.func.2_internal,
                               n=n, beta_estimate=beta_estimate, 
                               cali_prob = cali_prob, param_est = param_est,
                               data_case = data_case, data_ctrl = data_ctrl,
                               index.b = index.b))
    
    Numerator_1   <- numerator.1-numerator.2*numerator.3
    
    Numerator_2.1 <- B.case%*%(exp(t(beta_mle)%*%B.case)*
                                 sum(sapply(1:(num_match+1), function(x)exp(t(beta_mle)%*%B.cbind[,x]))))
    
    Numerator_2.2 <- c(exp(t(beta_mle)%*% B.case))*
      matrix(rowSums(sapply(1:(num_match+1), 
                            function(x)B.cbind[,x]*c(exp(t(beta_mle)%*% B.cbind[,x])))), ncol = 1)
    
    Denominator_2 <-  (sum(sapply(1:(num_match+1), 
                                  function(x)exp(t(beta_mle)%*% B.cbind[,x]))))^2
    
    Numerator_2   <- (Numerator_2.1/Denominator_2- Numerator_2.2/Denominator_2)/Denominator_2
    
    Numerator_all <- c(Numerator_1)*Numerator_2
    
    final.value <- (Numerator_all/L_sj)/L_sj
    
    x_b[, index.b] <- final.value
  }
  
  return( x_b )
}

#### returns for all strata (in list form) for one study
psi_beta_x_b_sj         <- function(beta_estimate, cali_prob, param_est,
                                    data_case, data_ctrl, #ncase = 1000, num_match = 1,
                                    calib_method){
  ncase <- length(unique(cali_prob[,"strat"]))
  beta_x_b_sj <- list()
  P           <- ncol(cali_prob) - 1 # calib_prob contains probs for each category + stratum number
  
  for(n in 1:ncase){
    cali_prob.temp <- cali_prob[cali_prob[, "strat"] == n, ]
    num_match      <- nrow(cali_prob.temp) - 1
    p_eval <- expand.grid(replicate( num_match + 1, 0:(P-1), simplify=FALSE))
    #### sapply returns a vector by column instead of a matrix 
    if(calib_method=="Full"){
      beta_x_b_sj.temp <- rowSums(sapply(split(p_eval, seq(nrow(p_eval))), psi_beta_x_b_sj_pi_full,
                                         n=n, beta_estimate=beta_estimate, 
                                         cali_prob = cali_prob, param_est = param_est,
                                         data_case = data_case, data_ctrl = data_ctrl))
    } else if(calib_method=="Internalized"){
      beta_x_b_sj.temp <- rowSums(sapply(split(p_eval, seq(nrow(p_eval))), psi_beta_x_b_sj_pi_internal,
                                         n=n, beta_estimate=beta_estimate, 
                                         cali_prob = cali_prob, param_est = param_est,
                                         data_case = data_case, data_ctrl = data_ctrl))
    }
    
    beta_x_b_sj[[n]] <- matrix(beta_x_b_sj.temp, nrow = length(beta_estimate), 
                               ncol = nrow(param_est), byrow = FALSE)
  }
  return(beta_x_b_sj)
}
mini.func.3             <- function(p, n, beta_estimate, cali_prob,
                                    data_case, data_ctrl, x.star){
  cali_prob.temp <- cali_prob[cali_prob[, "strat"] == n, ]
  num_match      <- nrow(cali_prob.temp) - 1
  p        <- unlist(p)
  p.case   <- p[1]  
  p.ctrl   <- p[-1] 
  
  #### beta_mle are the estimated parameters for beta_X 
  beta_mle <- matrix(beta_estimate, nrow = length(beta_estimate),
                     ncol = 1)
  #### B matrix for case
  B.case   <- matrix(0, nrow = nrow(beta_mle), ncol = 1)
  if(p.case!=0){
    B.case[p.case,] <- 1
  }
  B.ctrl   <- matrix(0, nrow = nrow(beta_mle), ncol = num_match)
  for(i in 1:num_match){
    if(p.ctrl[i]!=0){
      B.ctrl[p.ctrl[i],i] <- 1
    }
  }
  B.cbind  <- cbind(B.case, B.ctrl)
  
  temp_case <- data_case[data_case[,"strat"] == n,]
  temp_ctrl <- data_ctrl[data_ctrl[,"strat"] == n,]
  temp_cali_prob <- cali_prob.temp
  
  prt.1    <- prod(sapply(1:(num_match+1), 
                          function(x)temp_cali_prob[x, p[x]+1])) 
  useful_q <- sum(sapply(1:(num_match+1), function(x)exp(t(beta_mle)%*% B.cbind[,x])))
  prt.2    <- B.case[x.star,]*c(exp(t(beta_mle)%*%B.case))*useful_q
  prt.3_1  <- sum((sapply(1:(num_match+1), 
                          function(x)B.cbind[x.star,x]*c(exp(t(beta_mle)%*% B.cbind[,x])))))
  prt.3_2  <- c(exp(t(beta_mle)%*%B.case))
  prt.3    <- prt.3_1*prt.3_2
  
  final    <- prt.1*( (prt.2/c(useful_q)) - (prt.3/c(useful_q)) )/c(useful_q)
  return(final)
}
psi_beta_x_beta_x_sj_pi <- function(p, n, beta_estimate, cali_prob,
                                    data_case, data_ctrl){
  #### n corrsponds to the n-th case in the study
  ####    (or the n-th stratum)
  #### p is the p_2, \ldots, p_{M_{j(s)}+1} indicator, with
  #### first element be the case and the remaining one be the ctrls
  cali_prob.temp <- cali_prob[cali_prob[, "strat"] == n, ]
  num_match      <- nrow(cali_prob.temp) - 1
  P              <- ncol(cali_prob) - 1 # calib_prob contains probs for each category + stratum number
  
  p_eval <- expand.grid(replicate( num_match + 1, 0:(P-1), simplify=FALSE))
  L_sj   <- sum(sapply(split(p_eval, seq(nrow(p_eval))), lik_nk_fun_strata,
                       n = n, beta_mle0 = beta_estimate,
                       cali_prob = cali_prob))
  
  #### notice category starts from 0
  p        <- unlist(p)
  p.case   <- p[1]  
  p.ctrl   <- p[-1] 
  #### beta_mle are the estimated parameters for beta_X 
  beta_mle <- matrix(beta_estimate, nrow = length(beta_estimate),
                     ncol = 1)
  #### B matrix for case
  B.case   <- matrix(0, nrow = nrow(beta_mle), ncol = 1)
  if(p.case!=0){
    B.case[p.case,] <- 1
  }
  B.ctrl   <- matrix(0, nrow = nrow(beta_mle), ncol = num_match)
  for(i in 1:num_match){
    if(p.ctrl[i]!=0){
      B.ctrl[p.ctrl[i],i] <- 1
    }
  }
  B.cbind  <- cbind(B.case, B.ctrl)
  
  temp_case <- data_case[data_case[,"strat"] == n,]
  temp_ctrl <- data_ctrl[data_ctrl[,"strat"] == n,]
  temp_cali_prob <- cali_prob.temp
  W.case    <- temp_case$W
  W.ctrl    <- temp_ctrl$W
  W.combined <- c(W.case, W.ctrl)
  #### object to store the second order derivative
  x_x      <- matrix(NA, nrow = length(beta_estimate),
                     ncol = length(beta_estimate))
  useful_q <- sum(sapply(1:(num_match+1), function(x)exp(t(beta_mle)%*% B.cbind[,x])))
  for(x.prime in 1:length(beta_estimate)){
    for(x.star in 1:length(beta_estimate)){
      prt.1    <- prod(sapply(1:(num_match+1), 
                              function(x)temp_cali_prob[x, p[x]+1])) 
      
      prt.2    <- B.case[x.prime,1]*B.case[x.star,1]*c(exp(t(beta_mle)%*%B.case))*useful_q
      
      prt.3_1  <- B.case[x.prime,1]*c(exp(t(beta_mle)%*%B.case))
      prt.3_2  <- sum(sapply(1:(num_match+1), 
                             function(x)B.cbind[x.star,x]*exp(t(beta_mle)%*%B.cbind[,x]) ))
      prt.3    <- prt.3_1*prt.3_2
      
      prt.4_1  <- B.case[x.star,1]*c(exp(t(beta_mle)%*%B.case))
      prt.4_2  <- sum(sapply(1:(num_match+1), 
                             function(x)B.cbind[x.prime,x]*exp(t(beta_mle)%*%B.cbind[,x]) ))
      prt.4    <- prt.4_1*prt.4_2
      
      prt.5_1  <- sum(sapply(1:(num_match+1), 
                             function(x)B.cbind[x.prime,x]*B.cbind[x.star,x]*c(exp(t(beta_mle)%*%B.cbind[,x]))))
      prt.5_2  <- c(exp(t(beta_mle)%*%B.case))
      prt.5    <- prt.5_1*prt.5_2
      
      numerator.1 <- prt.1*(prt.2 + prt.3 - prt.4 - prt.5)*L_sj*(useful_q^2)
      
      prt.6    <- prt.1
      prt.7    <- B.case[x.prime,]*c(exp(t(beta_mle)%*%B.case))*useful_q
      prt.8_1  <- sum((sapply(1:(num_match+1), 
                              function(x)B.cbind[x.prime,x]*c(exp(t(beta_mle)%*%B.cbind[,x])))))
      prt.8_2  <- c(exp(t(beta_mle)%*%B.case))
      prt.8    <- prt.8_1*prt.8_2
      numerator.2_1 <- prt.6*(prt.7 - prt.8)
      
      prt.9    <- 2*useful_q
      prt.10   <- sum(sapply(1:(num_match+1), 
                             function(x)B.cbind[x.star,x]*exp(t(beta_mle)%*%B.cbind[,x])))*L_sj
      numerator.2_2 <- prt.9 * prt.10
      
      prt.11   <- useful_q^2
      prt.12   <- sum(sapply( split(p_eval, seq(nrow(p_eval))), mini.func.3,
                              n=n, beta_estimate=beta_estimate, 
                              cali_prob = cali_prob,
                              data_case = data_case, data_ctrl = data_ctrl,
                              x.star = x.star))
      numerator.2_3 <- prt.11*prt.12
      
      numerator.2   <- numerator.2_1*(numerator.2_2+numerator.2_3)
      
      Numerator     <- numerator.1 - numerator.2
      
      Denominator   <- (L_sj^2)*(useful_q^4)
      x_x[x.prime,x.star]  <- Numerator/Denominator
    }
  }
  return(x_x)
}
psi_beta_x_beta_x_sj    <- function(beta_estimate, cali_prob,
                                    data_case, data_ctrl){
  ncase <- length(unique(cali_prob[,"strat"]))
  betax_betax_sj <- list()
  P              <- ncol(cali_prob) - 1 # calib_prob contains probs for each category + stratum number
  
  for(n in 1:ncase){
    cali_prob.temp <- cali_prob[cali_prob[, "strat"] == n, ]
    num_match      <- nrow(cali_prob.temp) - 1
    p_eval <- expand.grid(replicate( num_match + 1, 0:(P-1), simplify=FALSE))
    #### sapply returns a vector by column instead of a matrix 
    betax_betax.temp <- rowSums(sapply(split(p_eval, seq(nrow(p_eval))), psi_beta_x_beta_x_sj_pi,
                                       n=n, beta_estimate=beta_estimate, 
                                       cali_prob = cali_prob, 
                                       data_case = data_case, data_ctrl = data_ctrl))
    betax_betax_sj[[n]] <- matrix(betax_betax.temp, nrow = length(beta_estimate), 
                                  ncol = length(beta_estimate), byrow = FALSE)
  }
  return(betax_betax_sj)
}

##################################################################
######################  Main Function  ###########################
##################################################################
ME_categorical          <- function(calib_method, P, mydata, S_set){
  # Parameters@calib_method: Either "Full" for full calibration method
  ##                         or "Internalized" for Internalized calibration method
  # Parameters@P:     Number of categories for the biomarker
  # Parameters@S_set: Number of studies for the pooling project
  # Input:            mydata - a list containing the data
  ## mydata[[s]]: contains the data for the s-th study
  
  ## mydata[[s]]$one_lab_case:          Data for the cases in the s-th study
  #### mydata[[s]]$one_lab_case$Y:      Binary disease outcome
  #### mydata[[s]]$one_lab_case$W:      Continuous local laboratory measurements
  #### mydata[[s]]$one_lab_case$strat:  stratum indicator, indicating which stratum this participant belongs to
  #### mydata[[s]]$one_lab_case$id:     id identifier
  ##### Remark: strat and id in the calibration subset need to be sorted
  ##### Both variables start with 1 and increase with step size 1
  
  ### mydata[[s]]$one_lab_ctrl:         Data for the matched controls in the s-th study
  #### mydata[[s]]$one_lab_ctrl$Y:      Binary disease outcome
  #### mydata[[s]]$one_lab_ctrl$W:      Continuous local laboratory measurements
  #### mydata[[s]]$one_lab_ctrl$strat:  Stratum indicator, indicating which stratum this participant belongs to
  #### mydata[[s]]$one_lab_ctrl$id:     id identifier
  ##### Remark: strat and id in the calibration subset need to be sorted
  ##### @strat should be matched to the cases in mydata[[s]]$one_lab_case
  ##### @id should start with nrow(mydata[[s]]$one_lab_case)+1 and increase with step size 1
  
  ### mydata[[s]]$one_reference:         Calibration subset
  #### mydata[[s]]$one_reference$Y:      Binary disease outcome
  #### mydata[[s]]$one_reference$W:      Continuous local laboratory measurements
  #### mydata[[s]]$one_reference$X_cate: Categorical reference laboratory measurements 
  #### mydata[[s]]$one_reference$strat:  stratum indicator, indicating which stratum this participant belongs to
  #### mydata[[s]]$one_reference$id:     id identifier
  ##### Remark: strat and id in the calibration subset do not need to be sorted
  ##### The categorical biomarker measurement X_cate should be recorded in the format of 0, 1, 2, 3,...
  
  # Output: Estimate_coef: estimated regression parameters in the main model
  #         SE: Standard errors of the estimated regression coefficients
  for(s in 1:S_set){
    one_reference <- mydata[[s]]$one_reference
    one_lab_case  <- mydata[[s]]$one_lab_case
    one_lab_ctrl  <- mydata[[s]]$one_lab_ctrl
    #### Fit the calibration model
    cali_model <- multinom( X_cate ~ W, data = one_reference, maxit = 100000)
    calib_estimate <- as.matrix(summary(cali_model)$coefficients)
    if(calib_method=="Full"){
      #### Full calibration: predict everyone
      cali_prob     <- predict(cali_model, 
                               newdata = data.frame( W = c(one_lab_case$W,one_lab_ctrl$W)), 
                               type = 'probs')
      cali_prob      <- cbind(cali_prob, strat = c(one_lab_case$strat, one_lab_ctrl$strat))
      cali_prob_ref  <- predict(cali_model, newdata = one_reference, type='probs')
    }else{
      cali_prob      <- predict(cali_model, 
                                newdata = data.frame( W = c(one_lab_case$W,one_lab_ctrl$W)), 
                                type = 'probs')
      cali_prob      <- cbind(cali_prob, strat = c(one_lab_case$strat, one_lab_ctrl$strat))
      cali_prob      <- cbind(cali_prob, id = c(one_lab_case$id, one_lab_ctrl$id))
      #internal validation study
      cali_prob[cali_prob[,"id"]%in%one_reference$id,1:P] <- 0
      for(iii in 0:(P-1)){
        id.c <- one_reference[one_reference$X_cate==iii, "id"]
        cali_prob[ cali_prob[,"id"] %in% id.c, iii+1] <- 1
      }
      cali_prob_ref <- cali_prob[cali_prob[,"id"] %in% one_reference$id, 1:(P+1)]
      cali_prob     <- cali_prob[,1:(P+1)] # remove id
    }
    mydata[[s]]   <- list(cali_prob      = cali_prob, 
                          cali_prob_ref  = cali_prob_ref,
                          one_reference  = one_reference, 
                          one_lab_case   = one_lab_case ,
                          one_lab_ctrl   = one_lab_ctrl, 
                          calib_estimate = calib_estimate)
  }
  
  cali_prob_lab <- lapply(mydata,"[[", 1) 
  beta_estimate <- optim(rep(0,P-1),lik_fun,
                         cali_prob_lab=cali_prob_lab,
                         method="L-BFGS-B",lower=-2,upper=2)
  # beta_estimate <- optim(c(0,0),lik_fun,
  #                        cali_prob_lab=cali_prob_lab,ncase=ncase_set,
  #                        method="L-BFGS-B",lower=-2,upper=2)
  beta_estimate <- beta_estimate$par
  #### Four parameters for the calibration model in each study
  A <- matrix(0, nrow = 2*(P-1)*S_set + length(beta_estimate),
              ncol = 2*(P-1)*S_set + length(beta_estimate))
  
  B <- matrix(0, nrow = 2*(P-1)*S_set + length(beta_estimate),
              ncol = 2*(P-1)*S_set + length(beta_estimate))
  for(s in 1:S_set){
    mydata_temp <- mydata[[s]]
    reference   <- mydata_temp$one_reference
    
    param_est   <- mydata_temp$calib_estimate
    cali_prob   <- mydata_temp$cali_prob
    data_case   <- mydata_temp$one_lab_case
    
    data_ctrl   <- mydata_temp$one_lab_ctrl
    
    case_ctrl   <- rbind(data_case, data_ctrl)
    case_ctrl$C <- case_ctrl$id%in%reference$id
    #### In case may sample both from cases and ctrls
    data_case$C <- data_case$id%in%reference$id
    data_ctrl$C <- data_ctrl$id%in%reference$id
    
    one_lab_main   = case_ctrl
    reference$strat <- NA
    for(i in reference$id){
      reference[reference$id==i,"strat"] <- one_lab_main[one_lab_main$id==i, "strat"]
    }
    
    ### Fill up A matrix 
    first_a <- psi_a_q_sj(reference = reference, param_est = param_est, 
                          cali_prob = cali_prob)
    first_b <- psi_b_q_sj(reference = reference, param_est = param_est, 
                          cali_prob = cali_prob)
    first_x <- psi_beta_x_sj(beta_estimate = beta_estimate, 
                             cali_prob = cali_prob)
    first_bind <- rbind(first_a, first_b, first_x)
    #### mark which rows in A should be filled by this study
    index.ina  <- c((2*(P-1)*(s-1)+1):(2*(P-1)*s), (2*(P-1)*(S_set)+1):nrow(A))
    count1 <- 1
    for(i in index.ina){
      count2 <- 1
      for(j in index.ina){
        temp <- sum(first_bind[count1,]*first_bind[count2,])
        A[i,j] <- A[i,j] + temp
        count2 <- count2 + 1
      }
      count1 <- count1 + 1
    }
    ### Fill up B matrix 
    
    ### The parts only involving a and b
    second_a  <- psi_double_a_q_s(reference = reference, 
                                  param_est = param_est)
    second_b  <- psi_double_b_q_s(reference = reference, 
                                  param_est = param_est)
    second_ab <- psi_a_b_q_s(reference = reference, 
                             param_est = param_est)
    index.ina_cali  <- c((2*(P-1)*(s-1)+1):(2*(P-1)*s))
    a.half <- index.ina_cali[1:(length(index.ina_cali)/2)]
    b.half <- index.ina_cali[(length(index.ina_cali)/2+1):(length(index.ina_cali))]
    B[a.half, a.half] <- B[a.half, a.half] + second_a
    B[b.half, b.half] <- B[b.half, b.half] + second_b
    B[a.half, b.half] <- B[a.half, b.half] + second_ab
    B[b.half, a.half] <- B[b.half, a.half] + second_ab
    
    #### The parts involving psi_x and derivative for a b and x (The last two rows)
    second_x_a <- psi_beta_x_a_sj(beta_estimate = beta_estimate, 
                                  cali_prob = cali_prob, param_est = param_est,
                                  data_case = data_case, data_ctrl = data_ctrl, 
                                  calib_method = calib_method)
    second_x_b <- psi_beta_x_b_sj(beta_estimate = beta_estimate, 
                                  cali_prob = cali_prob, param_est = param_est,
                                  data_case = data_case, data_ctrl = data_ctrl,
                                  calib_method = calib_method)
    second_x_x <- psi_beta_x_beta_x_sj(beta_estimate = beta_estimate, 
                                       cali_prob = cali_prob,
                                       data_case = data_case, data_ctrl = data_ctrl)
    
    ### which columns should be filled in B for one study
    ncase <- length(unique(cali_prob[,"strat"]))
    rows.for.x   <- (2*(P-1)*(S_set)+1):nrow(B)
    B[rows.for.x, a.half]    <- matrix(rowSums(sapply(1:ncase, function(x)second_x_a[[x]] )), 
                                       nrow = length(rows.for.x), ncol = length(a.half), byrow = FALSE) + B[rows.for.x, a.half]
    B[rows.for.x, b.half]    <- matrix(rowSums(sapply(1:ncase, function(x)second_x_b[[x]] )), 
                                       nrow = length(rows.for.x), ncol = length(b.half), byrow = FALSE) + B[rows.for.x, b.half]
    B[rows.for.x,rows.for.x] <- matrix(rowSums(sapply(1:ncase, function(x)second_x_x[[x]] )), 
                                       nrow = length(rows.for.x), ncol = length(rows.for.x), byrow = FALSE) + B[rows.for.x,rows.for.x]
  }
  sandwich.var <- try(solve(B)%*%A%*%t(solve(B)))
  sd.coef <- sqrt(diag(sandwich.var)[ (2*(P-1)*(S_set)+1):nrow(B)])
  return(list(Estimate_coef = beta_estimate,
              SE = sd.coef))
}

#### Working example
#### prepare data
S_set          <- 3
sd_beta0       <- 0.1
alpha_all      <- list()
alpha_all[[1]] <- matrix(c(-13.79246,  0.3258804,  
                           -29.82904,  0.6323898), byrow=T,nrow=2)
alpha_all[[2]] <- matrix(c(-13.51673,  0.3203838,   
                           -29.66921,  0.6348230),byrow=T,nrow=2)
alpha_all[[3]] <- matrix(c(-13.52193,  0.3246280,   
                           -29.70576,  0.6426687),byrow=T,nrow=2)
sd_w_all       <- c(16, 16, 16)
mu_w_all       <- c(33.86566, 41.32044, 47.76033)


beta_mle       <- (-log(1.5)/2)
num_match      <- 1
ncase          <- 200
ncase_set      <- ncase
nval           <- 50
n_v            <- nval

mydata <- list()
for(s in 1:S_set){
  #### Generate categorical X directly from W
  lab_data<-function(which_s, ncase, num_match, 
                     alpha_all, sd_w_all, beta, sd_beta0,
                     mu_w_all){
    #parameters
    #set.seed(seed.seed)
    beta  <- beta 
    #nctrl <- ncase*num_match
    alpha <- alpha_all[[which_s]]
    sd_w  <- sd_w_all[which_s]
    mu_w  <- mu_w_all[[which_s]]
    repeat{
      case <- data.frame(matrix(vector(), 0, 3, 
                                dimnames=list(c(), c("Y","W","X_cate"))),
                         stringsAsFactors=F)
      ctrl <- data.frame(matrix(vector(), 0, 3, 
                                dimnames=list(c(), c("Y","W","X_cate"))),
                         stringsAsFactors=F)
      for(stratum in 1:ncase){
        triples.temp <- data.frame(matrix(vector(), 0, 3, 
                                          dimnames=list(c(), c("Y","W","X_cate"))),
                                   stringsAsFactors=F)
        nr <- 1
        beta0_sj <- rnorm(1, -1, sd_beta0)
        while(sum(triples.temp$Y) < 1 | sum(1-triples.temp$Y) < num_match){
          #error-prone
          ##  Q1: local lab measurements are generated with same dist across labs
          W <- rnorm(1, mu_w, sd_w)
          
          ##  Q2: Calibration parameters all the same across labs
          #true
          prob_X <- c(1/(1+exp(alpha[1,1]+alpha[1,2]*W)+exp(alpha[2,1]+alpha[2,2]*W)),
                      exp(alpha[1,1]+alpha[1,2]*W)/(1+exp(alpha[1,1]+alpha[1,2]*W)+exp(alpha[2,1]+alpha[2,2]*W)),
                      exp(alpha[2,1]+alpha[2,2]*W)/(1+exp(alpha[1,1]+alpha[1,2]*W)+exp(alpha[2,1]+alpha[2,2]*W)))
          
          X_cate <- sample(c(0,1,2), 1, replace = FALSE, prob = prob_X)
          
          #### Q3:
          #### \beta_{0,sj} are all set to be 1, instead of beta_{0,sj}\sim N(0,sigma^2)
          
          prob_Y <- exp(beta0_sj + X_cate*beta)/(1 + exp(beta0_sj + X_cate*beta))
          Y      <- sample(c(1,0),1,prob=c(prob_Y,1-prob_Y))
          triples.temp[nr,] <- c(Y, W, X_cate)
          nr     <- nr + 1
        }
        case.temp <- triples.temp[which(triples.temp$Y==1)[1],]
        ctrl.temp <- triples.temp[which(triples.temp$Y==0)[1:num_match],]
        
        case       <- rbind(case, case.temp)
        ctrl       <- rbind(ctrl, ctrl.temp)
      }
      case$strat <- 1:ncase
      case$id    <- 1:ncase
      ctrl$strat <- rep(1:ncase, each = num_match)
      ctrl$id    <- (ncase + 1):(ncase*(num_match + 1))
      ### reorder ctrls
      index.ctrl <- sapply(1:num_match, function(x)seq(from = x, to = ncase*num_match, by = num_match))
      index.ctrl <- c(index.ctrl)
      ctrl       <- ctrl[index.ctrl,]
      if (sum(ctrl$X_cate==0) > 0 & sum(ctrl$X_cate==1) > 0 & sum(ctrl$X_cate==2) > 0){
        break
      }
    }
    return(list(case,ctrl))
  }
  one_lab <- lab_data(which_s = s, ncase = ncase, num_match = num_match,
                      alpha_all = alpha_all, sd_w_all = sd_w_all, mu_w_all = mu_w_all,
                      beta = beta_mle, sd_beta0 = sd_beta0)
  
  one_lab_case            <- one_lab[[1]]
  row.names(one_lab_case) <- c(1:nrow(one_lab_case))
  one_lab_ctrl            <- one_lab[[2]]
  row.names(one_lab_ctrl) <- c((nrow(one_lab_case)+1):(nrow(one_lab_case)+nrow(one_lab_ctrl)))
  #### Sample calibration subset based on quantiles of W
  one_lab_ctrl$qunatil    <- with(one_lab_ctrl,cut(one_lab_ctrl$W, quantile(one_lab_ctrl$W, prob = seq(0, 1, by=0.1)),
                                                   labels=factor(1:10), include.lowest=TRUE))
  ### Q4: When determining the calibration subset, 
  ###      ctrls are sampled by quantile of W
  repeat {
    one_reference <- NULL
    for (k in 1:10){
      one_reference <- rbind(one_lab_ctrl[sample(which(one_lab_ctrl$qunatil==k),n_v/10),], one_reference)
    }
    if (sum(one_reference$X_cate == 0)>0 & sum(one_reference$X_cate == 1)>0 & sum(one_reference$X_cate == 2)>0){
      break
    }
  }
  one_lab_ctrl <- one_lab_ctrl[,-6]
  data.temp <- list(one_reference  = one_reference, 
                    one_lab_case   = one_lab_case ,
                    one_lab_ctrl   = one_lab_ctrl)
  mydata[[s]] <- data.temp
}

#### Implementing the Full calibration method
full <- ME_categorical(P=3, mydata=mydata, S_set=3, calib_method = "Full")
# returns the point estimates
full$Estimate_coef
# returns the standard errors for the estimated coefficients
full$SE

#### Implementing the Internalized calibration method
internalized <- ME_categorical(P=3, mydata=mydata, S_set=3, calib_method = "Internalized")
internalized$Estimate_coef
internalized$SE
