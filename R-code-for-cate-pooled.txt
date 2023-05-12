####################################################################################
#                                                                                  #
#                                R Code                                            #
#                                                                                  #
#      Statistical Methods for Analysis of Combined Categorical Biomarker Data     #
#                            from Multiple Studies                                 #
#                                                                                  #
#                   For any question please contact Chao Cheng                     #                                                                                #    
#                         Email: cqplus@126.com                                    #
#                                                                                  #
#                                                                                  #
####################################################################################

#################################################################################################
#
# There are two sections, the 1st section is the functions, the 2nd is an illustrative example
#
#################################################################################################

#############################################################################################
#   
# SECTION 1: FUNCTIONS
#
#############################################################################################

##########################################################################
# R packages
##########################################################################
library(fGarch)


##########################################################################
##########################################################################
# calculate the gradient
##########################################################################
# INPUT
# f: traget function
# x0: point where the gradient is to build
# heps: step size
##########################################################################
# RETURN
# the numerical gradient of target function
##########################################################################
mygrad=function (f, x0,heps = 1e-6, ...) {
  if (!is.numeric(x0)) 
    stop("Argument 'x0' must be a numeric value.")
  fun <- match.fun(f)
  f <- function(x) fun(x, ...)
  p =length(f(x0))
  n <- length(x0)
  hh <- rep(0, n)
  gr <- matrix(0,nrow=n,ncol=p)
  for (i in 1:n) {
    hh[i] <- heps
    gr[i,] <- (f(x0 + hh) - f(x0 - hh))/(2 * heps)
    hh[i] <- 0
  }
  return(gr)
}

##########################################################################
##########################################################################
# logit function: f(x) = e^x/(1+e^x)
##########################################################################
# INPUT
# x: point where we need to calculate its logit value
##########################################################################
# RETURN
# the logit value of x
##########################################################################
logit=function(x) {exp(x)/(1+exp(x))}



##############################################################################################
##############################################################################################
# generate simulated data
# P(Y=1)=logit(alpha0 + beta1*X.cate + beta2*W),
# where X.cate is the categorical biomarker data generated from W based on some cut-off points
# beta2=0 denotes that W doesn't affect Y directly
##############################################################################################
# INPUT
# OR:         odd ratios (>0)
# prevalence: disease prevalence (0%~50%)
# beta2:      the coefficient of W
# cut:        cut-off points
# type:       1 for RSCS design and 2 for COCS design
# errorType:  the distrbution of the errors in X's model, "uniform" for uniform distribution,
#             "normal" for normal distribution, "skew" for skew normal distribution
# n1:         Number of individuals for each studies
##############################################################################################
# RETURN The simualted disease-biomarker data.
#        $data     the simulated data
#        $sigma    real sigma
#        $beta     real beta
##############################################################################################
data.sim=function(OR=0.5,prevalence=0.05,beta2=0,cut=5,type=1,errorType="normal",n1=1000){
  N = 15000
  M=5
  var_x = 300
  var_l = runif(M+1,150,250)
  alpha0_l = rnorm(M,60,sqrt(5))
  alpha1_l = 2
  
  mat=matrix(0,ncol = 6,nrow = M*N)
  G = length(cut)
  cate.mat = matrix(0,ncol=G,nrow=M*N)
  colnames(mat) = c("W","X","HC","HL","group","X.cut")
  for (i in (1:M)) {
    alpha0.now=alpha0_l[i]
    W   = rnorm(N,mean=0,sd=sqrt(100))
    if (errorType=="normal") {
      X   = rnorm(N,mean = alpha0.now + alpha1_l*W ,sd=sqrt(var_x))
    } else if(errorType=="uniform") {
      X   = alpha0.now+alpha1_l*W + runif(N,min=-sqrt(3*var_x),max=sqrt(3*var_x))
    } else {
      X   = alpha0.now+alpha1_l*W + fGarch::rsnorm(N, mean = 0, sd = sqrt(var_x), xi = 3)
    }
    H_C = rnorm(N,mean=X,sd=sqrt(var_l[1]))
    H_L = rnorm(N,mean=X,sd=sqrt(var_l[i+1]))
    mat[((i-1)*N+1):(i*N),1] = W
    mat[((i-1)*N+1):(i*N),2] = X
    mat[((i-1)*N+1):(i*N),3] = H_C
    mat[((i-1)*N+1):(i*N),4] = H_L
    mat[((i-1)*N+1):(i*N),5] = i
    mat[((i-1)*N+1):(i*N),6] = apply(X> rep(1,N) %*% t(cut),1,sum)
    cate.mat[((i-1)*N+1):(i*N),] = ((X> rep(1,N) %*% t(cut)) - (X>rep(1,N) %*% t(c(cut[-1],Inf))))
  }
  cut.ex=c(-Inf,cut[-length(cut)])
  par.coe=mu.beta=matrix(0,ncol=G+1,nrow=M)
  for (i in (1:M)) {
    par.coe[i,] = c(log(OR),beta2); 
    mu.beta[i,] = c(pnorm(c(cut[-1],+Inf),mean=alpha0_l[i] + alpha1_l*0,sd=sqrt(alpha1_l*alpha1_l*5+var_x))-pnorm(cut,mean=alpha0_l[i] + alpha1_l*0,sd=sqrt(alpha1_l*alpha1_l*5+var_x)),0) 
  }
  beta0=rep(0,M)
  prob.vector=rep(0,M*N)
  for (i in (1:M)) {
    beta0[i] = log(prevalence/(1-prevalence)) - sum(par.coe[i,]*mu.beta[i,])
    prob.vector[((i-1)*N+1):(i*N)] = logit(
      beta0[i] + t(log(OR)) %*% t(cate.mat[((i-1)*N+1):(i*N),]) + beta2*mat[((i-1)*N+1):(i*N),"W"])
  }
  y <- as.numeric(runif(M*N) < prob.vector)
  mat=cbind(mat,cate.mat)
  if (type==1) {
    n.cali= as.integer(n1/10)
    data=matrix(0,ncol = dim(mat)[2]+1,nrow = M*n1)
    is.calibration=rep(0,M*n1)
    for (i in (1:M)) {
      y.now   = y[((i-1)*N+1):(i*N)]
      mat.now = mat[((i-1)*N+1):(i*N),]
      loc = sample(1:N,n1,replace=F)
      data[((i-1)*n1+1):(i*n1),-1] = mat.now[loc,]
      data[((i-1)*n1+1):(i*n1),1] = y.now[loc]
      a = sample(1:n1,n.cali,replace=F)
      is.calibration[(i-1)*n1+a]=1
    }
  } else if (type==2){
    data=matrix(0,ncol = dim(mat)[2]+1,nrow = M*n1)
    is.calibration=rep(0,M*n1)
    n.cali= as.integer(n1/10)
    for (i in (1:M)) {
      y.now   = y[((i-1)*N+1):(i*N)]
      mat.now = mat[((i-1)*N+1):(i*N),]
      loc = sample(1:N,n1,replace=F)
      data[((i-1)*n1+1):(i*n1),-1] = mat.now[loc,]
      data[((i-1)*n1+1):(i*n1),1] = y.now[loc]
      a = sample(which(data[((i-1)*n1+1):(i*n1),1]==0),n.cali,replace=F)
      is.calibration[(i-1)*n1+a]=1
    }
  } else {
    data=matrix(0,ncol = dim(mat)[2]+1,nrow = M*1000)
    is.calibration=rep(0,M*1000)
    n.cali= as.integer(n1/10)
    for (i in (1:M)) {
      y.now   = y[((i-1)*N+1):(i*N)]
      mat.now = mat[((i-1)*N+1):(i*N),]
      loc1 = sample(which(y.now==1),min(sum(y.now),500),replace=F)
      if ( (1000-length(loc1)) > length( which(y.now==0))) {
        loc2 = sample(which(y.now==0),1000-length(loc1),replace=T)
      } else {
        loc2 = sample(which(y.now==0),1000-length(loc1),replace=F)
      }
      data[((i-1)*1000+1):(i*1000),-1] = mat.now[c(loc1,loc2),]
      data[((i-1)*1000+1):(i*1000),1] = y.now[c(loc1,loc2)]
      a = sample(which(data[((i-1)*1000+1):(i*1000),1]==0),100,replace=F)
      is.calibration[(i-1)*1000+a]=1
    }
  }
  data=cbind(data,is.calibration)
  colnames(data)=c("Y","W","X","HC","HL","group","X.cut",paste("cate",1:G,sep="_"),"calibration")
  data[which(data[,"calibration"]==0),"HC"]=NA
  list(data=data,sigma=c(var_x,var_l),beta=c(beta0,log(OR),beta2))
}

##############################################################################################
##############################################################################################
# calculate beta0 in data.sim function
##############################################################################################
# RETURN: beta0
##############################################################################################
beta0.f=function(prevalence=0.05,cut=0,param=c(1,1),mu=rep(0,2),sigma=diag(2),type=1) {
  ex = pnorm(cut,mean=mu[2],sd=sqrt(sigma[2,2]))
  log(prevalence/(1-prevalence)) - param[2]*ex - param[1]*mu[1]
}


##############################################################################################
##############################################################################################
# calculate beta, as well as the first stage in our two-stage iteration method
##############################################################################################
##############################################################################################
# INPUT
# data: the disease-biomarker data frame
# Wname: the name of additional covariates in X's model
# sigma: sigma 
##############################################################################################
# RETURN: theta
############################################################################################## 
get.theta=function(data=mydata,Wname=c("W"),sigma=rep(1,7)) {
  n.W = length(Wname)
  M=length(sigma)-2
  var_x=sigma[1];var_l=sigma[-1]
  calibration = data[which(data[,"calibration"]==1),]
  ordinary    = data[which(data[,"calibration"]==0),]
  ## DH
  DH=array(0,dim=c(2,M+n.W,dim(calibration)[1]))
  DM=array(0,dim=c(1,M+n.W,dim(ordinary)[1]))
  for (i in (1:M)) {
    DH[1,i,]   <- calibration[,"group"]==i
    DH[2,i,]   <- calibration[,"group"]==i
    DM[1,i,]   <- ordinary[,"group"]==i
  }
  if (is.null(Wname)==F) {
    for (i in (1:n.W)) {
      DH[1,M+i,] = DH[2,M+i,] = calibration[,Wname[i]]
      DM[1,M+i,] = ordinary[,Wname[i]]
    }
  }
  
  ## sigma
  sigma1=array(var_x,dim=c(2,2,dim(calibration)[1]))
  sigma2=array(var_x,dim=c(1,1,dim(ordinary)[1]))
  sigma1[1,1,] = sigma1[1,1,] + var_l[1]
  for (i in (1:M)) {
    sigma1[2,2,which(calibration[,"group"]==i)] = var_x + var_l[i+1]
    sigma2[1,1,which(ordinary[,"group"]==i)] = var_x + var_l[i+1]
  }
  HC = array(0,dim=c(2,1,dim(calibration)[1]))
  HM = array(ordinary[,"HL"],dim=c(1,1,dim(ordinary)[1]))
  HC[1,1,]=calibration[,"HC"];HC[2,1,]=calibration[,"HL"]
  HC0 = HC[,1,]
  ## transfer 3 demensional array to list
  DH = lapply(seq(dim(DH)[3]), function(x) DH[ , , x])
  DM = lapply(seq(dim(DM)[3]), function(x) DM[ , , x])
  sigma1 = lapply(seq(dim(sigma1)[3]), function(x) sigma1[ , , x])
  sigma2_I = as.list(1/sigma2)
  HC = lapply(seq(dim(HC)[3]), function(x) HC[ , , x])
  HM=as.list(HM)
  A = apply(mapply(function(a,b,c) t(a) %*% solve(b) %*% c,a=DH,b=sigma1,c=HC),1,sum)
  B = apply(mapply(function(a,b,c) a * b * c,a=DM,b=sigma2_I,c=HM),1,sum)
  
  C = apply(mapply(function(a,b,c) t(a) %*% solve(b) %*% c,a=DH,b=sigma1,c=DH),1,sum)
  C = matrix(C,ncol=M+n.W,nrow=M+n.W)
  D = apply(mapply(function(a,b,c) b * a %*% t(c),a=DM,b=sigma2_I,c=DM),1,sum)
  D=matrix(D,ncol=M+n.W,nrow=M+n.W)
  ## theta
  theta = as.vector(solve(C+D) %*% (A+B))
  if (is.null(Wname)) {
    names(theta) = c(paste("alpha_0",1:M,sep=""))
  } else {
    names(theta) = c(paste("alpha_0",1:M,sep=""),paste("tau",1:n.W,sep=""))
  }
  
  ## ec & em
  EC = HC0 - mapply(function(x) x%*%theta ,x=DH)
  EM = ordinary[,"HL"] - mapply(function(x) x%*%theta,x=DM)
  list(theta=theta,ec=EC,em=EM,mats=list(A,B,C,D))
}

##############################################################################################
##############################################################################################
# calculate sigma, as well as the second stage in our two-stage iteration method
##############################################################################################
# INPUT
# data: the disease-biomarker data frame
# theta.out: the output of get.theta function
##############################################################################################
# RETURN: sigma
##############################################################################################
get.sigma=function(data,theta.out) {
  theta = theta.out$theta; EC=theta.out$ec;EM=theta.out$em
  M=max(data[,"group"]) 
  calibration = data[which(data[,"calibration"]==1),]
  ordinary    = data[which(data[,"calibration"]==0),]
  ## define A and B
  A=B=list()
  E=matrix(0,nrow=3,ncol=dim(calibration)[1])
  E[1,]=EC[1,]^2;E[2,]=EC[1,]*EC[2,];E[3,]=EC[2,]^2
  Y1 = matrix(0,nrow=M+2,ncol=dim(calibration)[1])
  Y2 = matrix(0,nrow=M+2,ncol=dim(ordinary)[1])
  X1 = X2 = matrix(0,ncol=M+2,nrow=M+2)
  for (i in (1:M)) {
    loc1 = which(calibration[,"group"]==i)
    loc2 = which(ordinary[,"group"]==i)
    
    A[[i]] = matrix(0,ncol=3,nrow=M+2)
    A[[i]][1,]=1;A[[i]][2,]=c(1,0,0)
    A[[i]][i+2,3]=1
    B[[i]] = A[[i]][,3] 
    
    Y1[,loc1] = A[[i]] %*% E[,loc1]
    Y2[,loc2] = B[[i]] %*% t(EM[loc2]^2)
    X1[,1]   = X1[,1] + length(loc1) * A[[i]] %*% rep(1,3)
    X1[,2]   = X1[,2] + length(loc1) * A[[i]] %*% c(1,0,0)
    X1[,i+2] = length(loc1)*A[[i]] %*% c(0,0,1)
    X2   = X2 + length(loc2) * crossprod(t(B[[i]]))
  }
  Y = apply(Y1,1,sum,na.rm=T) + apply(Y2,1,sum,na.rm=T)
  output = as.vector(solve(X1+X2) %*% Y)
  names(output)=c("sigma_x",paste("sigma_",0:M,sep=""))
  list(output,mats=list(A,B))
}

##############################################################################################
##############################################################################################
# calculate beta, as well as the second stage in our two-stage iteration method
##############################################################################################
# RETURN: 
# estimated beta by naive, cut-off calibration and exact calibration method
# $par.n  estimated beta by naive method
# $par.c  estimated beta by cut-off calibration method
# $par.p  estimated beta by exact calibration method
# $data   other data generated by this function
##############################################################################################
get.beta=function(data=mydata,theta,sigma,cut=5,includeW=T,Wname=c("W"),Zname="w2") {
  
  n.cate = length(cut)+1
  M=length(sigma)-2
  sigma=abs(sigma)
  n.W = length(Wname)
  
  if (is.null(Zname)) {
    includeZ=FALSE
    par0=rep(1,M+n.cate-1)
  } else {
    includeZ=TRUE
    Zname=Zname
    n.Z=length(Zname)
    par0=rep(1,M+n.cate-1+n.Z)
  }
  
  tau = theta[-(1:(M))]
  cut.ex = c(-Inf,cut,Inf)
  N = dim(data)[1]
  
  w1 = sigma[1]/(sigma[1] + sigma[-(1:2)])
  w2 = sigma[1]/(sigma[1]+sigma[2]*sigma[-c(1,2)]/(sigma[2]+sigma[-c(1,2)]))
  w2.in = sigma[-c(1:2)]/(sigma[2]+sigma[-c(1,2)])
  
  alpha = theta[1:M]
  sigma.hat.X1 = sigma[1]*sigma[-c(1,2)]/(sigma[1] + sigma[-c(1,2)])
  sigma.hat.X2 = w2 * sigma[2]*sigma[-c(1,2)]/(sigma[2]+sigma[-c(1,2)])
  
  Xe = H.mean = rep(0,dim(data)[1])
  prob.vec = matrix(0,ncol=length(cut)+1,nrow=dim(data)[1])
  for (i in (1:M)) {
    loc1 = which(data[,"group"]==i & data[,"calibration"]==0)
    loc2 = which(data[,"group"]==i & data[,"calibration"]==1)
    Xe[loc1] = w1[i]*(data[loc1,"HL"]) + (1-w1[i])*( alpha[i] + data[loc1,Wname] %*% matrix(tau,ncol=1))
    Xe[loc2] = w2[i]*(  w2.in[i]*(data[loc2,"HC"]) + (1-w2.in[i])*(data[loc2,"HL"])  ) + (1-w2[i])*(alpha[i] + data[loc2,Wname] %*% matrix(tau,ncol=1))
    H.mean[loc1] = data[loc1,"HL"]
    H.mean[loc2] = (data[loc2,"HL"] + data[loc2,"HC"])/2
    for (j in (1:(length(cut)+1))) {
      prob.vec[loc1,j] = pnorm(cut.ex[j+1],mean=Xe[loc1],sd = sqrt(sigma.hat.X1[i]))-pnorm(cut.ex[j],mean=Xe[loc1],sd = sqrt(sigma.hat.X1[i]))
      prob.vec[loc2,j] = pnorm(cut.ex[j+1],mean=Xe[loc2],sd = sqrt(sigma.hat.X2[i]))-pnorm(cut.ex[j],mean=Xe[loc2],sd = sqrt(sigma.hat.X2[i]))
    }
  }
  Xe.cut <- ((Xe> rep(1,N) %*% t(cut)) - (Xe>rep(1,N) %*% t(c(cut[-1],Inf))))
  colnames(prob.vec) = paste("prob_",1:n.cate,sep="")
  colnames(Xe.cut)=paste("Xe.cut",1:length(cut),sep="_")
  
  HL.cut <- ((H.mean> rep(1,N) %*% t(cut)) - (H.mean>rep(1,N) %*% t(c(cut[-1],Inf))))
  colnames(HL.cut)=paste("HL.cut",1:length(cut),sep="_")
  
  data = cbind(data,Xe,Xe.cut,HL.cut,prob.vec)
  
  loc.b1=(M+1):(M+length(cut))
  if (includeZ==T) {
    loc.b2 = (M+length(cut)+1):(M+length(cut)+length(Zname))
  }
  ### calibration method
  loglike=function(par) {
    l=0
    for (i in (1:M)) {
      loc = which(data[,"group"]==i)
      if (includeZ==T) {
        A= par[i] + Xe.cut[loc,] %*% matrix(par[loc.b1],ncol=1) + data[loc,Zname] %*% matrix(par[loc.b2],ncol=1)
      }else {
        A= par[i] + Xe.cut[loc,] %*% matrix(par[loc.b1],ncol=1)
      }
      l=l+sum( log( exp(data[loc,"Y"]*A)/(1+exp(A))) ) 
    }
    -l
  }
  ### naive method
  loglike0=function(par) {
    l=0
    for (i in (1:M)) {
      loc = which(data[,"group"]==i)
      if (includeZ==T) {
        A= par[i] + HL.cut[loc,] %*% matrix(par[loc.b1],ncol=1) + data[loc,Zname] %*% matrix(par[loc.b2],ncol=1)
      }else {
        A= par[i] + HL.cut[loc,] %*% matrix(par[loc.b1],ncol=1)
      }
      l=l+sum( log( exp(data[loc,"Y"]*A)/(1+exp(A))) ) 
    }
    -l
  }
  n.cate=length(cut)+1
  ### probabilitic method
  loglike1=function(par) {
    beta.extend = c(0,par[loc.b1])
    par.int=rep(0,N)
    for (i in (1:M)) {
      loc = which(data[,"group"]==i)
      par.int[loc] = par[i]
    }
    if (includeZ) {
      A = par.int %*% t(rep(1,n.cate)) + rep(1,N) %*% t(beta.extend) + data[,Zname] %*% matrix(par[loc.b2],ncol=1) %*% t(rep(1,n.cate))
    } else {
      A = par.int %*% t(rep(1,n.cate)) + rep(1,N) %*% t(beta.extend)
    }
    z = apply(prob.vec * exp((data[,"Y"] %*% t(rep(1,n.cate)))*(A))/(1+exp(A)),1,sum)
    -sum(log(z),na.rm=T)
  }
  dlogit1=function(par) {
    beta.extend = c(0,par[loc.b1])
    par.int=rep(0,N)
    for (i in (1:M)) {
      loc = which(data[,"group"]==i)
      par.int[loc] = par[i]
    }
    if (includeZ) {
      A= par.int %*% t(rep(1,n.cate)) + rep(1,N) %*% t(beta.extend) + data[,Zname] %*% matrix(par[loc.b2],ncol=1) %*% t(rep(1,n.cate))
    } else {
      A = par.int %*% t(rep(1,n.cate)) + rep(1,N) %*% t(beta.extend)
    }
    Y.mat = data[,"Y"] %*% t(rep(1,n.cate))
    ilogit = 1/apply(((exp(Y.mat*A)/(1+exp(A)))*prob.vec),1,sum)
    dl=matrix(0,ncol=length(par),nrow=N)
    for (i in (1:M)) {
      loc = which(data[,"group"]==i)
      Y.mat.now = Y.mat[loc,];A.now=A[loc,]
      D = (  Y.mat.now * exp(Y.mat.now * A.now) + 
               (Y.mat.now-1)* exp((Y.mat.now+1)*A.now)   )  /  (1+exp(A.now))^2
      dl[loc,i] = apply(D * prob.vec[loc,],1,sum) * ilogit[loc]
    }
    dl[,(M+1):(M+n.cate-1)]=prob.vec[,-1]*(Y.mat[,-1] * exp(Y.mat[,-1]*A[,-1])
                                           + (Y.mat[,-1]-1)*exp((Y.mat[,-1]+1)*A[,-1])) /  ((1+exp(A[,-1]))^2) * (ilogit %*% matrix(1,ncol=n.cate-1,nrow=1))
    if (includeZ==T) {
      for (i in (1:length(Zname))) {
        W.mat = data[,Zname[i]] %*% t(rep(1,n.cate)) 
        G = ((W.mat*Y.mat * exp(Y.mat * A) + W.mat*(Y.mat-1)*exp((Y.mat+1)*A)) / (1+exp(A))^2)
        dl[,M+n.cate+i-1] = apply((G*prob.vec),1,sum)*ilogit
      }
    }
    -apply(dl,2,sum,na.rm=T)
  }
  par.n = optim(par=par0,fn=loglike0,method="BFGS")$par
  par.c = optim(par=par.n,fn=loglike,method="BFGS")$par
  par.p = optim(par=par.c,fn=loglike1,gr=dlogit1,method="BFGS",control = list(reltol=1e-16))$par
  
  if(includeZ==T) {
    names(par.c) = names(par.n) = names(par.p)= c(paste("beta_0",1:M,sep="")
                                                  ,paste("beta_x",1:length(cut),sep=""),paste("beta",Zname,sep="_"))
  } else {
    names(par.c) = names(par.n) = names(par.p)= c(paste("beta_0",1:M,sep="")
                                                  ,paste("beta_x",1:length(cut),sep=""))
  }
  list(par.c=par.c,par.n=par.n,par.p=par.p,data=data)
}

##############################################################################################
##############################################################################################
# calculate 95% confidence interval of \hat(\beta)
##############################################################################################
# RETURN: the 95% confidence interval of \hat(\beta)
#         $beta1 for the esimated confidence interval of beta1
#         $sd.n the estimated standard errors of naive method
#         $sd.c the estimated standard errors of cut-off calibartion method
#         $sd.p the estimated standard errors of exact calibration method
#         $sd.n1 the estimated standard errors of naive method with hessian matrix
#         $sd.c1 the estimated standard errors of cut-off method with hessian matrix
#         $sd.p1 the estimated standard errors of exact method with hessian matrix
##############################################################################################
get.ci=function(data=out.beta,theta=theta.mat[3,],sigma=sigma.mat[3,],includeW=F
                ,Wname="W",Zname="w2",cut=cut,ErrorMat=sigma.now$ErrorMat) {

  mat = data$data;par.n=data$par.n;par.c=data$par.c;par.p=data$par.p
  M = sum(substr(names(par.n),6,6)=="0")
  N = dim(mat)[1]
  n.cate = sum(substr(names(par.n),6,6)=="x")+1
  n.beta = length(beta);n.theta=length(theta);n.sigma=length(sigma)
  prob.vec = mat[,which(substr(colnames(mat),1,4)=="prob")]
  HL.cut = mat[,which(substr(colnames(mat),1,6)=="HL.cut")]
  cut.ex=c(-Inf,cut,Inf)
  n.W = length(Wname)
  if (is.null(Zname)) {
    includeZ=FALSE
  } else {
    includeZ=TRUE
    Zname=Zname
    n.Z=length(Zname)
  }
  
  ##### the theta
  calibration = mat[which(mat[,"calibration"]==1),]
  ordinary    = mat[which(mat[,"calibration"]==0),]
  ## DH
  DH=array(0,dim=c(2,M+n.W,dim(calibration)[1]))
  DM=array(0,dim=c(1,M+n.W,dim(ordinary)[1]))
  for (i in (1:M)) {
    DH[1,i,]   <- calibration[,"group"]==i
    DH[2,i,]   <- calibration[,"group"]==i
    DM[1,i,]   <- ordinary[,"group"]==i
  }
  if (is.null(Wname)==F) {
    for (j in (1:M)) {
      loc1 = calibration[,"group"]==j
      loc2 = ordinary[,"group"]==j
      for (i in (1:n.W)) {
        DH[1,M+i,] = DH[2,M+i,] = calibration[,Wname[i]]*loc1
        DM[1,M+i,] = ordinary[,Wname[i]]*loc2
      }
    }
  } 
  HC = array(0,dim=c(2,1,dim(calibration)[1]))
  HM = array(ordinary[,"HL"],dim=c(1,1,dim(ordinary)[1]))
  HC[1,1,]=calibration[,"HC"];HC[2,1,]=calibration[,"HL"]
  HC0 = HC[,1,]
  ## transfer 3 demensional array to list
  DH = lapply(seq(dim(DH)[3]), function(x) DH[ , , x])
  DM = lapply(seq(dim(DM)[3]), function(x) DM[ , , x])
  HC = lapply(seq(dim(HC)[3]), function(x) HC[ , , x])
  HM=as.list(HM)
  d.l1=function(sigma,theta,type=1) {
    tau = theta[-c(1:(M))]
    ErrorMat = rep(0,3*M)
    for (i in (1:M)) {
      ErrorMat[(i-1)*3+1]=  sigma[1] + sigma[2]
      ErrorMat[(i-1)*3+2]=  sigma[1] 
      ErrorMat[(i-1)*3+3]=  sigma[1] + sigma[i+2]
    }
    
    ## sigma
    sigma1=array(0,dim=c(2,2,dim(calibration)[1]))
    sigma2=array(0,dim=c(1,1,dim(ordinary)[1]))
    for (i in (1:M)) {
      loc1=which(calibration[,"group"]==i)
      loc2=which(ordinary[,"group"]==i)
      sigma1[1,1,loc1] = ErrorMat[3*(i-1)+1]
      sigma1[1,2,loc1] = sigma1[2,1,loc1] = ErrorMat[3*(i-1)+2]
      sigma1[2,2,loc1] = ErrorMat[3*(i-1)+3]
      
      sigma2[1,1,loc2] = ErrorMat[3*(i-1)+3]
    }
    sigma1 = lapply(seq(dim(sigma1)[3]), function(x) sigma1[ , , x])
    sigma2_I = as.list(1/sigma2)
    
    A = mapply(function(a,b,c) t(a) %*% solve(b) %*% (c-a%*%theta),a=DH,b=sigma1,c=HC)
    B = mapply(function(a,b,c) a * b * (c-a%*%theta),a=DM,b=sigma2_I,c=HM)
    C = matrix(0,ncol=N,nrow=length(theta))
    
    l1=l2=0
    for (i in (1:M)) {
      loc1 = which(mat[,"calibration"]==1 & mat[,"group"]==i)
      loc2 = which(mat[,"calibration"]==0 & mat[,"group"]==i)
      C[,loc1] = A[,(l1+1):(l1+length(loc1))]
      C[,loc2] = B[,(l2+1):(l2+length(loc2))]
      l1 = (l1+length(loc1))
      l2 = (l2+length(loc2))
    }
    if (type==1) {
      D = mapply(function(a,b) t(a) %*% solve(b) %*% a,a=DH,b=sigma1)
      E = mapply(function(a,b) b*a %*% t(a),a=DM,b=sigma2_I)
      d2 = matrix(apply(D,1,sum),nrow=length(theta))+matrix(apply(E,1,sum),nrow=length(theta))
      list(Ut = t(C), U1t = d2)
    } else {
      U1s = apply(t(C),2,sum)
    }
  }
  
  XX = d.l1(sigma=sigma,theta=theta,type=1)
  Ut = XX$Ut; Q1t = XX$U1t
  Q1s = t(mygrad(d.l1,x0=sigma,type=2,theta=theta,0.000001))
  
  Q1 = cbind(Q1t,Q1s,matrix(0,nrow=dim(Q1t)[1],ncol=length(par.n)))
  
  
  ####### sigma
  A_sigma = matrix(0,ncol=N,nrow=3*length(sigma)); B_sigma = matrix(0,ncol=N,nrow=length(sigma))
  A2_sigma = B2_sigma = matrix(0,ncol=N,nrow=length(sigma)^2)
  for (i in (1:M)) {
    up = c(1,1,rep(0,M))
    mid = c(1,rep(0,M+1))
    down = rep(0,M+2)
    down[1] = 1; down[2+i] = 1
    loc1 = which(calibration[,"group"]==i)
    loc2 = which(ordinary[,"group"]==i)
    
    A_sigma[,loc1]= c(up,mid,down)
    B_sigma[,loc2]= down
    A2_sigma[,loc1] = as.vector(matrix(c(up,mid,down),ncol=3) %*% t(matrix(c(up,mid,down),ncol=3)))
    B2_sigma[,loc2] = down %*% t(down)
  }
  Q2s = matrix(apply(A2_sigma,1,sum) + apply(B2_sigma,1,sum),ncol=M+2)
  
  d.l2=function(theta,sigma,type=1) {
    tau = theta[-c(1:(M))]
    
    EC = mapply(function(a,b) b-a%*%theta,a=DH,b=HC)
    EM = mapply(function(a,b) b-a%*%theta,a=DM,b=HM)
    E1 = rbind(EC[1,]^2,EC[1,]*EC[2,],EC[2,]^2)
    E2 = EM^2
    d2 = matrix(0,nrow=length(sigma),ncol=N)
    for (i in (1:M)) {
      loc1 = which(calibration[,"group"]==i)
      loc2 = which(ordinary[,"group"]==i)
      sig1 = c( sigma[1] + sigma[2],sigma[1], sigma[1]+sigma[2+i])
      sig2 = sig1[3]
      
      E1.now = E1[,loc1] - sig1 %*% t(rep(1,length(loc1)))
      E2.now = E2[loc2] - rep(sig2,length(loc2))
      
      up = c(1,1,rep(0,M))
      mid = c(1,rep(0,M+1))
      down = rep(0,M+2)
      down[1] = 1; down[2+i] = 1
      
      A1 = matrix(c(up,mid,down),ncol=3) %*% E1.now 
      A2 = matrix(down,ncol=1) %*% t(E2.now)
      
      d2[,which(mat[,"calibration"]==1 & mat[,"group"]==i)] = A1
      d2[,which(mat[,"calibration"]==0 & mat[,"group"]==i)] = A2
    }
    if (type==1) {
      t(d2)
    } else {
      apply(d2,1,sum)
    }
  }
  Us=d.l2(theta=theta,sigma=sigma,type=1)
  Q2t= t(mygrad(d.l2,x0=theta,type=2,sigma=sigma,0.000001))
  Q2 = cbind(Q2t,Q2s,matrix(0,nrow=length(sigma),ncol=length(par.n)))
  
  ## naive or cut-off calibration method
  d.loglike=function(par,type=1,cut.name="HL.cut") {
    theta = par[1:n.theta];sigma=par[-(1:n.theta)][1:n.sigma];beta = par[-(1:(n.theta+n.sigma))]
    alpha=theta[1:M]
    tau = theta[-(1:(M))]
    
    sigma=abs(sigma)
    beta2 = beta[-c(1:(M+n.cate-1))]
    if (cut.name=="HL.cut") {
      if (is.null(dim(HL.cut))) {
        Xe.cut = matrix(HL.cut,ncol=1)
      } else {
        Xe.cut = HL.cut
      }
    } else {
      w1 = 1 * sigma[1] /(1 * sigma[1] + sigma[-c(1,2)])
      sigma.hat.X1 = w1 * sigma[-c(1,2)]
      w2 = sigma[1]/( sigma[1] + sigma[2]*sigma[-c(1:2)]/(1 * sigma[-c(1:2)] + 1 * sigma[2]) )                
      w2.in = 1 * sigma[-c(1:2)] /( 1 * sigma[-c(1:2)] + 1 * sigma[2])
      sigma.hat.X2 = w2 * sigma[2]*sigma[-c(1:2)]/(1 * sigma[-c(1:2)] + 1 * sigma[2])
      Xe = rep(0,N)
      for (i in (1:M)) {
        loc1 = which(mat[,"group"]==i & mat[,"calibration"]==0)
        loc2 = which(mat[,"group"]==i & mat[,"calibration"]==1)
        Xe[loc1] = w1[i]*(mat[loc1,"HL"] ) + (1-w1[i])*( alpha[i] + mat[loc1,Wname] %*% matrix(tau,ncol=1))
        Xe[loc2] = w2[i]*(  w2.in[i]*(mat[loc2,"HC"]) + (1-w2.in[i])*(mat[loc2,"HL"])  ) + (1-w2[i])*(alpha[i] + mat[loc2,Wname] %*% matrix(tau,ncol=1))
      }
      Xe.cut <- ((Xe> rep(1,N) %*% t(cut)) - (Xe>rep(1,N) %*% t(c(cut[-1],Inf))))
      colnames(Xe.cut)=paste("Xe.cut",1:length(cut),sep="_")
    }
    Logit = rep(0,N)
    if(includeZ) {
      psi = matrix(0,ncol=M+n.cate-1+length(Zname),nrow=N)
    } else {
      psi = matrix(0,ncol=M+n.cate-1,nrow=N)
    }
    for (i in (1:M)) {
      loc = which(mat[,"group"]==i)
      xe = Xe.cut[loc,]
      w  = mat[loc,Zname]
      y  = mat[loc,"Y"]
      if (includeZ) {
        Logit[loc] = y - logit(xe %*% matrix(beta[(M+1):(M+n.cate-1)],ncol=1) + w %*% matrix(beta2,ncol=1) + beta[i])
      } else {
        Logit[loc] = y - logit(xe %*% matrix(beta[(M+1):(M+n.cate-1)],ncol=1) + beta[i])
      }
      psi[loc,i] = Logit[loc]
    }
    psi[,(M+1):(M+n.cate-1)] = (Logit %*% matrix(1,ncol=n.cate-1,nrow=1)) * Xe.cut
    if (includeZ) {
      psi[,(M+n.cate):length(beta)] = (Logit %*% matrix(1,ncol=length(Zname),nrow=1)) * mat[,Zname] 
    }
    if (type==1) {
      psi
    } else {
      apply(psi,2,sum)
    }
  }
  ####### naive method
  param.n=c(theta,sigma,par.n)
  Ub.n = d.loglike(param.n,type=1,cut.name = "HL.cut")
  Qb.n = mygrad(f=d.loglike,x0=param.n,
                type=2,0.000001,cut.name="HL.cut")
  ####### cut-off method
  param.c=c(theta,sigma,par.c)
  Ub.c = d.loglike(param.c,type=1,cut.name = "Xe.cut")
  Qb.c = mygrad(f=d.loglike,x0=param.n,
                type=2,0.000001,cut.name="Xe.cut")
  #### exact calibration method
  d.loglike2=function(par,type=1) {
    theta = par[1:n.theta];sigma=par[-(1:n.theta)][1:n.sigma];beta = par[-(1:(n.theta+n.sigma))]
    sigma=abs(sigma)
    alpha=theta[1:M]
    tau = theta[-(1:M)]
    
    w1 =  sigma[1] /( sigma[1] + sigma[-c(1,2)])
    sigma.hat.X1 = w1 * sigma[-c(1,2)]
    
    w2 = sigma[1]/( sigma[1] + sigma[2]*sigma[-c(1:2)]/( sigma[-c(1:2)] +  sigma[2]) )                
    w2.in = sigma[-c(1:2)] /(  sigma[-c(1:2)] +   sigma[2])
    sigma.hat.X2 = w2 * sigma[2]*sigma[-c(1:2)]/( sigma[-c(1:2)] +  sigma[2])
    
    Xe = rep(0,dim(mat)[1])
    prob.vec = matrix(0,ncol=length(cut)+1,nrow=dim(mat)[1])
    for (i in (1:M)) {
      loc1 = which(mat[,"group"]==i & mat[,"calibration"]==0)
      loc2 = which(mat[,"group"]==i & mat[,"calibration"]==1)
      Xe[loc1] = w1[i]*(mat[loc1,"HL"]) + (1-w1[i])*( alpha[i] + mat[loc1,Wname] %*% matrix(tau,ncol=1))
      Xe[loc2] = w2[i]*(  w2.in[i]*(mat[loc2,"HC"]) + (1-w2.in[i])*(mat[loc2,"HL"]) ) + (1-w2[i])*(alpha[i] + mat[loc2,Wname] %*% matrix(tau,ncol=1))
      for (j in (1:(length(cut)+1))) {
        prob.vec[loc1,j] = pnorm(cut.ex[j+1],mean=Xe[loc1],sd = sqrt(sigma.hat.X1[i]))-pnorm(cut.ex[j],mean=Xe[loc1],sd = sqrt(sigma.hat.X1[i]))
        prob.vec[loc2,j] = pnorm(cut.ex[j+1],mean=Xe[loc2],sd = sqrt(sigma.hat.X2[i]))-pnorm(cut.ex[j],mean=Xe[loc2],sd = sqrt(sigma.hat.X2[i]))
      }
    }
    colnames(prob.vec) = paste("prob_",1:n.cate,sep="")
    
    
    beta.extend = c(0,beta[(M+1):(M+n.cate-1)])
    par.int=rep(0,N)
    for (i in (1:M)) {
      loc = which(mat[,"group"]==i)
      par.int[loc] = beta[i]
    }
    if (includeZ) {
      A = par.int %*% t(rep(1,n.cate)) + rep(1,N) %*% t(beta.extend) + mat[,Zname] %*% matrix(beta[(M+n.cate):length(beta)],ncol=1) %*% t(rep(1,n.cate))
    } else {
      A = par.int %*% t(rep(1,n.cate)) + rep(1,N) %*% t(beta.extend)
    }
    
    Y.mat = mat[,"Y"] %*% t(rep(1,n.cate))
    ilogit = 1/apply(((exp(Y.mat*A)/(1+exp(A)))*prob.vec),1,sum)
    dl=matrix(0,ncol=length(beta),nrow=N)
    for (i in (1:M)) {
      loc = which(mat[,"group"]==i)
      Y.mat.now = Y.mat[loc,];A.now=A[loc,]
      D = (  Y.mat.now * exp(Y.mat.now * A.now) + 
               (Y.mat.now-1)* exp((Y.mat.now+1)*A.now)   )  /  (1+exp(A.now))^2
      dl[loc,i] = apply(D * prob.vec[loc,],1,sum) * ilogit[loc]
    }
    dl[,(M+1):(M+n.cate-1)]=prob.vec[,-1]*(Y.mat[,-1] * exp(Y.mat[,-1]*A[,-1])
                                           + (Y.mat[,-1]-1)*exp((Y.mat[,-1]+1)*A[,-1])   ) /  (1+exp(A[,-1]))^2 * (ilogit %*% matrix(1,ncol=n.cate-1,nrow=1))
    if (includeZ==T) {
      for (i in (1:length(Zname))) {
        W.mat = mat[,Zname[i]] %*% t(rep(1,n.cate)) 
        G = ((W.mat*Y.mat * exp(Y.mat * A) + W.mat*(Y.mat-1)*exp((Y.mat+1)*A)) / (1+exp(A))^2)
        dl[,M+n.cate+i-1] = apply((G*prob.vec),1,sum)*ilogit
      }
    }
    if (type==1) {
      dl
    } else {
      apply(dl,2,sum,na.rm=T)
    }
  }
  
  param.p=c(theta,sigma,par.p)
  Ub.p = d.loglike2(param.p,type=1)
  Qb.p = mygrad(f=d.loglike2,x0=param.p,
                type=2,0.000001)
  Qb.p[(n.theta+1):(n.theta+n.sigma),]=0
  ## covariance matrix
  ##### naive method
  Q.n = rbind(Q1,Q2,t(Qb.n))
  U.n.s = cbind(Ut,Us,Ub.n)
  U.n = t(U.n.s) %*% U.n.s
  cov.n = solve(Q.n) %*% U.n %*% t(solve(Q.n))
  
  Qb.n1 = Qb.n[-c(1:(n.sigma+n.theta)),]
  cov.n1 = -solve(Qb.n1)
  ## cut-off
  Q.c = rbind(Q1,Q2,t(Qb.c))
  U.c.s = cbind(Ut,Us,Ub.c)
  U.c = t(U.c.s) %*% U.c.s
  cov.c = solve(Q.c) %*% U.c %*% t(solve(Q.c))
  
  Qb.c1 = Qb.c[-c(1:(n.sigma+n.theta)),]
  cov.c1 = -solve(Qb.c1)
  ## exact method
  Q.p = rbind(Q1,Q2,t(Qb.p))
  U.p.s = cbind(Ut,Us,Ub.p)
  U.p = t(U.p.s) %*% U.p.s
  cov.p = solve(Q.p) %*% U.p %*% t(solve(Q.p))
  
  Qb.p0 = Qb.p
  Qb.p0[1:(n.theta+n.sigma),]=0
  Q.p0 = rbind(Q1,Q2,t(Qb.p0))
  U.p.s = cbind(Ut,Us,Ub.p)
  U.p = t(U.p.s) %*% U.p.s
  cov.p0 = solve(Q.p0) %*% U.p %*% t(solve(Q.p0))
  
  Q.p1 = Qb.p0[-c(1:(n.sigma+n.theta)),]
  cov.p1 = -solve(Q.p1)
  
  output=matrix(0,nrow=6,ncol=2*(n.cate-1))
  loc.beta = n.theta+n.sigma+M
  colName=NULL
  for (i in (1:(n.cate-1))) {
    r1 = sqrt(diag(cov.n))[loc.beta+i]
    r2 = sqrt(diag(cov.n1))[M+i]
    r3 = sqrt(diag(cov.c))[loc.beta+i]
    r4 = sqrt(diag(cov.c1))[M+i]
    r5 = sqrt(diag(cov.p))[loc.beta+i]
    r6 = sqrt(diag(cov.p1))[M+i]
    output[1,(2*i-1):(2*i)] = c(par.n[M+i]-1.96*r1,par.n[M+i]+1.96*r1)
    output[2,(2*i-1):(2*i)] = c(par.n[M+i]-1.96*r2,par.n[M+i]+1.96*r2)
    output[3,(2*i-1):(2*i)] = c(par.c[M+i]-1.96*r3,par.c[M+i]+1.96*r3)
    output[4,(2*i-1):(2*i)] = c(par.c[M+i]-1.96*r4,par.c[M+i]+1.96*r4)
    output[5,(2*i-1):(2*i)] = c(par.p[M+i]-1.96*r5,par.p[M+i]+1.96*r5)
    output[6,(2*i-1):(2*i)] = c(par.p[M+i]-1.96*r6,par.p[M+i]+1.96*r6)
    colName=c(colName,paste("beta_x",i,"_lower",sep=""),paste("beta_x",i,"_upper",sep=""))
  }
  colnames(output) = colName
  rownames(output) = c("CI-naive-sandwich","CI-naive-hessian","CI-cutoff-sandwich","CI-cutoff-hessian",
                       "CI-exact-sandwich","CI-exact-hessian")
  list(beta1=output,sd.n=sqrt(diag(cov.n)),sd.c=sqrt(diag(cov.c)),sd.p=sqrt(diag(cov.p)),
       sd.n1=sqrt(diag(cov.n1)),sd.c1=sqrt(diag(cov.c1)),sd.p1=sqrt(diag(cov.p1)))
}


##############################################################################################
##############################################################################################
# MAIN estimation function
#calculate theta, sigma, beta and corresponding S.E
##############################################################################################
# INPUT
# mat:               a data.frame or matrix containing all the information we need
# Y_name:            disease outcome
# HL_name:           local lab measurement
# HC_name:           central lab measurement
# group_name:        the label of the participating study (1,2,...,M)
# calibration_name:  whether selected into the calibration subset, 0 for no, 1 for yes
# W_name:            additional covariates in the X's model, NULL for none
# Z_name:            additional covariates in the logistic regression, NULL for none
# cut_points:        the cut-off points
# sandwich:          1 for sandwich method, 0 for hessian method
##############################################################################################
# RETURN: a list including estimates of theta, sigma, beta and conresponding standard errors
#         $theta       the estimates of theta
#         $sigma       the estimates of variance for the error terms in X's model
#         $beta.naive  the estimated beta by naive method
#         $beta.cutoff the estimated beta by cut-off calibration method
#         $beta.exact  the estimated beta by exact calibration method
#         $ci          the 95% confidence interval of log(OR) for above three methods
#         $se.naive    the S.E for all estimates by naive method
#         $se.cutoff   the S.E for all estimates by cut-off calibartion method
#         $se.exact    the S.E for all estimates by exact calibartion method
##############################################################################################

mainfunc=function(mat,Y_name,HL_name,HC_name,group_name,calibration_name,W_name,Z_name
                  ,cut_points,sandwich=TRUE) {
  options(warn=-1)
  M = max(mat[,group_name])
  data = cbind(Y=mat[,Y_name],HL=mat[,HL_name],HC=mat[,HC_name],group=mat[,group_name]
               ,calibration=mat[,calibration_name],mat[,W_name,drop=F],mat[,Z_name,drop=F])
  theta.mat=matrix(0,ncol=M+length(W_name),nrow=100);
  sigma.mat=matrix(0,ncol=M+2,nrow=100);
  if (is.null(W_name)) {
    Name_X = paste("alpha",1:M,sep="")
  } else {
    Name_X = c(paste("alpha",1:M,sep=""),W_name)
  }
  Name_sigma = c("sigma_x",  "sigma_0",  paste("sigma",1:M,sep="_"))
  i=1
  stop.cert=1
  while( stop.cert>0.01 & i<=99) {
    if (i==1) {
      out=get.theta(data=data,sigma=rep(1,M+2),Wname=W_name)
      theta.mat[i,]=out$theta
    } else {
      out=get.theta(data=data,sigma=sigma.now[[1]],Wname=W_name)
      theta.mat[i,]=out$theta
    }
    sigma.now=get.sigma(data=data,theta.out=out)
    sigma.mat[i,]=sigma.now[[1]]
    stop.cert=ifelse(i==1,1,sum(abs(theta.mat[i,]-theta.mat[i-1,]))+sum(abs(theta.mat[i,]-theta.mat[i-1,])))
    i=i+1
  }
  j=i-1
  includeW=ifelse(is.null(W_name),FALSE,TRUE)
  out.beta = get.beta(data=data,cut=cut_points,theta=theta.mat[j,],sigma=abs(sigma.mat[j,]),includeW=includeW,Wname=W_name,Zname=Z_name)
  ci = get.ci(data=out.beta,theta=theta.mat[j,],sigma=abs(sigma.mat[j,]),includeW=includeW,Wname=W_name,cut=cut_points,Zname=Z_name)
  
  res.theta = theta.mat[j,];names(res.theta) = Name_X
  res.sigma = sigma.mat[j,];names(res.sigma) = Name_sigma
  if (sandwich==1) {
    res.ci = ci$beta1[c(1,3,5),]
    res.sd.n = ci$sd.n
    names(res.sd.n) = c(Name_X,Name_sigma,names(out.beta$par.n))
    res.sd.c = ci$sd.c
    names(res.sd.c) = c(Name_X,Name_sigma,names(out.beta$par.c))
    res.sd.p = ci$sd.p
    names(res.sd.p) = c(Name_X,Name_sigma,names(out.beta$par.p))
  } else {
    res.ci = ci$beta1[c(2,4,6),]
    res.sd.n = c(ci$sd.n[1:(length(Name_X)+length(Name_sigma))],ci$sd.n1)
    names(res.sd.n) = c(Name_X,Name_sigma,names(out.beta$par.n))
    res.sd.c = c(ci$sd.n[1:(length(Name_X)+length(Name_sigma))],ci$sd.c1)
    names(res.sd.c) = c(Name_X,Name_sigma,names(out.beta$par.c))
    res.sd.p = c(ci$sd.n[1:(length(Name_X)+length(Name_sigma))],ci$sd.p1)
    names(res.sd.p) = c(Name_X,Name_sigma,names(out.beta$par.p))
  }
  list(theta=res.theta,sigma=res.sigma,beta.naive=out.beta$par.n,beta.cutoff=out.beta$par.c,beta.exact=out.beta$par.p,
       ci=res.ci,se.naive=res.sd.n,se.cutoff=res.sd.n,se.exact=res.sd.p)
}


#############################################################################################
#
# SECTION 2: AN ILLUSTRATIVE EXAMPLE
#
#############################################################################################

#############################################
#1, generate data
#############################################

### Odd ratios
OR = exp(c(0.5,1)*log(2))
### disease prevalence
prevalence=0.5
### Controls only calibration design
type = 2
### W does not affect Y directly; if W affect Y directly, you can set any value does not equal 0.
beta2=0
### cut-off points
cut_points = c(40,65)
### simulate data
data = data.sim(OR=OR,prevalence=prevalence,beta2=beta2,type=type,cut=cut_points,n1=1000)




#############################################
#2, parameter estimation
#############################################
results = mainfunc(mat = data$data                  # dataframe
                   ,Y_name = "Y"                    # disease outcome
                   ,HL_name = "HL"                  # local lab measurement
                   ,HC_name = "HC"                  # central lab measurement
                   ,group_name = "group"            # study label
                   ,calibration_name="calibration"  # calibration label
                   ,W_name = c("W")                 # covariates affecting X
                   ,Z_name = NULL                   # covariates affecting Y
                   ,cut_points=cut_points           # cut-off points
                   ,sandwich=0                      # 1 for sandwich method, 0 for hessian method
)

### theta
results$theta
### sigma
results$sigma
### beta by naive method and its conrsponding S.E
results$beta.naive
results$se.naive
### beta by cut-off method and its conrsponding S.E
results$beta.cutoff
results$se.cutoff
### beta by exact method and its conrsponding S.E
results$beta.exact
results$se.exact
### 95% confidence interval of beta_x
results$ci

