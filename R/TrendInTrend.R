# version 2020 Jan
#' An Odds Ratio Estimation Function
#' @description Estimate causal odds ratio (OR) given trends in exposure prevalence and outcome frequencies of stratified data. 
#' @param n11 A G by Tn matrix with n11[i,j] being the count of treated subjects with an event within group i at time j. The number of strata is G and the number of time intervals is Tn.
#' @param n10 A G by Tn matrix with n10[i,j] being the count of treated subjects without an event within group i at time j.
#' @param n01 A G by Tn matrix with n01[i,j] being the count of untreated subjects with an event within group i at time j.
#' @param n00 A G by Tn matrix with n00[i,j] being the count of untreated subjects without an event within group i at time j.
#' @param bnull Initial values for beta0, beta1, beta2 for the optimization algorithm. Default is (-10,0,0). It is suggested the initial value of beta0 be set as a small negative number (-4 or smaller) for the rare outcome model to be computationally stable.
#' @param n_explore Number of iterations in the optimization algorithm to stabilize the outputs. Default is 10.
#' @param noise_var The optimization algorithm is iterated n_explore times. Results from the previous iteration with added Gaussian noise are set as the starting values for the new iteration. Bigger noise_var indicates larger variance for the Gaussian noise, meaning more exploration during the iterations. Default is (1,1,0.5).
#' @param n_boot Number of bootstrap iterations to construct the confidence interval for the estimated odds ratio beta1. Default is 50.
#' @param alpha (1-alpha) is the significance level of the confidence interval. Default is 0.05.
#' @return \item{beta}{Maximum likelihood estimators (MLE) for beta0, beta1, beta2. Beta1 is the estimated treatment-event odds ratio. Because we conduct n_explore iterations, the set of parameters that is associated with the highest log likelihood is the output.} 
#' \item{CI_beta1}{1-alpha confidence interval for beta1.}
#' \item{ll}{Log likelihood evaluated at the MLE.}
#' \item{not_identified}{Equals 1 if the MLE is not identifiable or weakly identified. This could happen when there are multiple sets of parameters associated with the highest log likelihood, or the bootstrap confidence interval fails to cover the estimated beta1.}
#' @details This function estimates the odds ratio parameter beta1 in the subject-specific model in Ji et al. (2017)
#' \deqn{logit(E[Y(it)|Z(it), G(i), X(it)])=beta0+Z(it)*beta1+t*beta2+X(it)\gamma}
#'where \eqn{Z(it)} and \eqn{Y(it)} are the binary exposure and outcome variables for individual \eqn{i} at time \eqn{t}. 
#'There are three caveats regarding the implementation. First, the trend-in-trend design works better when there are substantial exposure trend differences across strata. If the exposure trend is roughly parallel across strata, the method may fail to converge. Second, we recommend running the OR function for multiple starting points to evaluate the stability of the optimization algorithm. Third, the bootstrap confidence interval may have slightly lower coverage probability than the nominal significance level 1-alpha.
#' @importFrom stats optim rnorm rbinom 
#' @importFrom stats glm predict quantile
#' @importFrom grDevices rainbow
#' @import nleqslv 
#' @import pracma 
#' @import pROC
#' @references Ji X, Small DS, Leonard CE, Hennessy S (2017). The Trend-in-trend Research Design for Causal Inference. Epidemiology 28(4), 529–536.  \cr
#' Ertefaie, A., Small, D., Ji, X., Leonard, C., Hennessy, S. (2018). Statistical Power for Trend-in-trend Design. Epidemiology 29(3), e21.\cr
#' Ertefaie, A., Small, D., Leonard, C., Ji, X., Hennessy, S. (2018). Assumptions Underlying the Trend-in-Trend Research Design. Epidemiology  29(6), e52-e53.
#' @examples
#' \donttest{
#' data <- GenData()
#' n11 <- data[[1]]
#' n10 <- data[[2]]
#' n01 <- data[[3]]
#' n00 <- data[[4]]
#' results <- OR(n11,n10,n01,n00)
#' }
#' @export
OR<-function(n11,n10,n01,n00,bnull=c(-10,0,0),n_explore=10,noise_var=c(1,1,0.5),n_boot=50,alpha=0.05){
  tmp<-OR_opt(n11,n10,n01,n00,bnull,noise_var,n_explore)
  para_index<-tmp$para_index
  not_identified<-tmp$not_identified
  h2<-tmp$h2
  ll<-tmp$ll
  if(not_identified==1){
    CI_boot<-NA
  }else{
    CI_boot<-var_boot(n11,n10,n01,n00,para_index,h2,bnull,noise_var,n_explore,n_boot,alpha)
    ind<-1*(para_index[2]> CI_boot[1] & para_index[2]<CI_boot[2]) # check whether CI covers beta1
    if(is.na(CI_boot[1]) | ind==0){
      warning("this odds ratio estimate may not be reliable.") 
      not_identified<-1
      CI_boot<-NA
    }else{
      CI_boot<-round(CI_boot,3)
    }
  }
  output<-list(beta=para_index[1:3],CI_beta1=CI_boot,ll=ll,not_identified=not_identified)
  output
}

OR_opt<-function(n11,n10,n01,n00,bnull,noise_var,n_explore){
  G<-dim(n11)[1]
  Tn<-dim(n11)[2]
  for(i in 1:G)
    for(j in 1:Tn)
      if(n11[i,j]+n10[i,j]==0)
        n10[i,j]<-1
  h2<-rep(NA,Tn)
  h2[1]<-0.001
  for(i in 2:Tn){
    h2[i]<-h2[i-1]*sum(n11[,i]+n10[,i])/sum(n11[,i]+n10[,i]+n01[,i]+n00[,i])*sum(n11[,i-1]+n10[,i-1]+n01[,i-1]+n00[,i-1])/sum(n11[,i-1]+n10[,i-1])
  }
  b0<-bnull
  G0<-rep(c(0.5,1,0.2),rep(G,3))
  start<-c(b0,G0)
  noise_var_tmp<-noise_var
  para.const<-list(n11,n10,n01,n00,h2,G,Tn)
  
  opt_res_est<-matrix(data=NA,nrow = n_explore,ncol = length(start))
  opt_res_var<-matrix(data=NA,nrow = n_explore,ncol = length(bnull))
  opt_res_ll<-matrix(data=NA,nrow=n_explore,ncol=1)
  for(i in 1:n_explore){
    mle1<-optim(start,fn=LL,gr=grr,para.const=list(n11,n10,n01,n00,h2,G,Tn),control=list("fnscale"=-1, maxit=10000), method = "Nelder-Mead",hessian=FALSE)
    converge<-mle1$convergence
    if(converge==0) start[1:length(bnull)]<-round(mle1$par[1:length(bnull)],3) 
    
    mle<-nleqslv(x=start,fn=grr,para.const=para.const,control=list(maxit=10000),method="Newton") # could also try Broyden
    if(abs(mle$x[2])>50) {start[1:length(bnull)]<-0} # stablizing step, if beta 1 is estimated too large, restart at 0
    else{
      start[1:length(bnull)]<-round(mle$x[1:3],3)
    }
    mle<-optim(start,fn=LL,gr=grr,para.const=list(n11,n10,n01,n00,h2,G,Tn),control=list("fnscale"=-1, maxit=10000), method = "Nelder-Mead",hessian=FALSE) 
    converge<-mle$convergence
    if(converge!=0){
      opt_res_est[i,]<-0
      opt_res_ll[i,1]<--Inf
      noise_var_tmp<-noise_var
    }
    if(converge==0){
      opt_res_est[i,]<-round(mle$par,3)
      opt_res_ll[i,1]<-round(LL(mle$par,para.const=list(n11,n10,n01,n00,h2,G,Tn)),2)
      noise_var_tmp<-noise_var_tmp/sqrt(i)
    }
    start[1:length(bnull)]<-rnorm(length(bnull),opt_res_est[i,],noise_var_tmp) # the initial value for the next iteration, which is last time outputs+noise
    #start[1:length(bnull)]<-rnorm(length(bnull),opt_res_est[i,],noise_var) # the initial value for the next iteration, which is last time outputs+noise
    }
  not_identified<-0
  if(length(unique(opt_res_ll[,1]))==1) {warning("Multiple solutions identified!")
    not_identified<-1
    para_index<-NA
    h2<-NA
    ll<-NA
  }else{
    index<-which(opt_res_ll[,1]==max(opt_res_ll[,1]))[1]
    para_index<-opt_res_est[index,]
    ll<-opt_res_ll[index,1]
  }
  output<-list(para_index=para_index,h2=h2,not_identified=not_identified,ll=ll)
  return(output)
}

var_boot<-function(n11,n10,n01,n00, para,h2,bnull,noise_var,n_explore,n_boot,alpha){
  G<-dim(n11)[1]
  Tn<-dim(n11)[2]
  beta0<-para[1]
  beta1<-para[2] # this is the estimated value from original
  beta2<-para[3]
  C1 <- para[4:(4+G-1)]
  C2 <- para[(4+G):(4+2*G-1)]
  C3 <- para[(4+2*G):(4+3*G-1)]
  n_Z1<-n11+n10
  n_Z0<-n01+n00
  res<-matrix(data=NA,nrow=n_boot,ncol=length(para))
  n11_b<-matrix(data=NA,nrow = G, ncol=Tn)
  n10_b<-matrix(data=NA,nrow = G, ncol=Tn)
  n01_b<-matrix(data=NA,nrow = G, ncol=Tn)
  n00_b<-matrix(data=NA,nrow = G, ncol=Tn)
  for(b in 1:n_boot){
    print(b)
    for(g in 1:G){
      for(t in 1:Tn){
        n11_b[g,t]<-rbinom(1,n_Z1[g,t],exp(beta0+beta1+beta2*t)*C1[g]/C3[g])
        n01_b[g,t]<-rbinom(1,n_Z0[g,t],exp(beta0+beta2*t)*(C2[g]-h2[t]*C1[g])/(1-h2[t]*C3[g]))
      }
    }
    n10_b<-n_Z1-n11_b
    n00_b<-n_Z0-n01_b
    tryCatch({
      tmp<-OR_opt(n11_b,n10_b,n01_b,n00_b,bnull,noise_var,n_explore)
      res[b,]<-tmp[[1]]},error=function(e) print(e))
  }
  CI<-quantile(res[,2],c(alpha/2,1-alpha/2),na.rm = TRUE)
  return(CI)
}

LL<-function(para,para.const){
  b<-rep(NA,3)
  b[1]<-para[1]
  b[2]<-para[2]
  b[3]<-para[3]
  n11<-para.const[[1]]
  n10<-para.const[[2]]
  n01<-para.const[[3]]
  n00<-para.const[[4]]
  h2<-para.const[[5]]
  G<-para.const[[6]]
  Tn<-para.const[[7]]
  C1 <- para[4:(4+G-1)]
  C2 <- para[(4+G):(4+2*G-1)]
  C3 <- para[(4+2*G):(4+3*G-1)]
  loglik<-0
  for(g in 1:G)
    for(t in 1:Tn){
      if(C1[g]<=0)
      {
        C1[g]<-exp(-10)
      }
      if(C3[g]<=0)
      {
        C3[g]<-exp(-10)
      }
      if((C2[g] - C1[g] * h2[t])<0)
      {
        C2[g]<-C1[g] * h2[t]+exp(-10)
      }
      if((1 - C3[g] * h2[t])<0)
      {
        C3[g]<-1/h2[t]-exp(-10)
      }
      if(exp(b[1]+b[3]*t)*(C2[g]-C1[g]*h2[t])/(1-C3[g]*h2[t])>1)
        b[1]<-log(1/(exp(b[3]*t)*(C2[g]-C1[g]*h2[t])/(1-C3[g]*h2[t])))-exp(-10)
      if(exp(b[1]+b[2]+b[3]*t)*C1[g]/C3[g]>1)
        b[1]<-log(1/(exp(b[2]+b[3]*t)*C1[g]/C3[g]))-exp(-10)
      
      loglik<-loglik+n11[g,t]*(b[1]+b[2]+b[3]*t) + n11[g,t]*(log(C1[g])-log(C3[g]))+
        n10[g,t]*log(1-exp(b[1]+b[2]+b[3]*t)*C1[g]/C3[g])+
        n01[g,t]*(b[1]+b[3]*t)+n01[g,t]*(log(C2[g]-C1[g]*h2[t])-log(1-C3[g]*h2[t]))+
        n00[g,t]*log(1-exp(b[1]+b[3]*t)*(C2[g]-C1[g]*h2[t])/(1-C3[g]*h2[t]))
    }
  loglik
}

grr<-function(para,para.const){
  b<-rep(NA,3)
  b[1]<-para[1]
  b[2]<-para[2]
  b[3]<-para[3]
  n11<-para.const[[1]]
  n10<-para.const[[2]]
  n01<-para.const[[3]]
  n00<-para.const[[4]]
  h2<-para.const[[5]]
  G<-para.const[[6]]
  Tn<-para.const[[7]]
  C1 <- para[4:(4+G-1)]
  C2 <- para[(4+G):(4+2*G-1)]
  C3 <- para[(4+2*G):(4+3*G-1)]
  gr_res<-numeric(length(para))
  for(g in 1:G){
    for(t in 1:Tn){
      if(C1[g]<=0)
      {
        C1[g]<-exp(-10)
      }
      if(C3[g]<=0)
      {
        C3[g]<-exp(-10)
      }
      if((C2[g] - C1[g] * h2[t])<0)
      {
        C2[g]<-C1[g] * h2[t]+exp(-10)
      }
      if((1 - C3[g] * h2[t])<0)
      {
        C3[g]<-1/h2[t]-exp(-10)
      }
      if(exp(b[1]+b[3]*t)*(C2[g]-C1[g]*h2[t])/(1-C3[g]*h2[t])>1)
        b[1]<-log(1/(exp(b[3]*t)*(C2[g]-C1[g]*h2[t])/(1-C3[g]*h2[t])))-exp(-10)
      if(exp(b[1]+b[2]+b[3]*t)*C1[g]/C3[g]>1)
        b[1]<-log(1/(exp(b[2]+b[3]*t)*C1[g]/C3[g]))-exp(-10)
      
      tmp1<-(C1[g]/C3[g]*exp(b[1]+b[2]+b[3]*t))
      tmp2<-exp(b[1]+b[3]*t)*(C2[g]-C1[g]*h2[t])/(1-C3[g]*h2[t])
      gr_res[1]<-gr_res[1]+n11[g,t]-n10[g,t]*tmp1/(1-tmp1)+n01[g,t]-n00[g,t]*tmp2/(1-tmp2)
      gr_res[2]<-gr_res[2]+n11[g,t]-n10[g,t]*tmp1/(1-tmp1)
      gr_res[3]<-gr_res[3]+t*(n11[g,t]-n10[g,t]*tmp1/(1-tmp1)+n01[g,t]-n00[g,t]*tmp2/(1-tmp2))
        
      gr_res[3+g]<-gr_res[3+g]+n11[g,t]/C1[g]-n10[g,t]*exp(b[1]+b[2]+b[3]*t)/C3[g]/(1-tmp1)-n01[g,t]*h2[t]/(C2[g]-h2[t]*C1[g])+
        n00[g,t]*exp(b[1]+t*b[3])*h2[t]/(1-h2[t]*C3[g])/(1-tmp2)
      gr_res[3+G+g]<-gr_res[3+G+g]+n01[g,t]/(C2[g]-h2[t]*C1[g])-n00[g,t]*exp(b[1]+b[3]*t)/(1-h2[t]*C3[g])/(1-tmp2)
      gr_res[3+2*G+g]<-gr_res[3+2*G+g]-n11[g,t]/C3[g]+n10[g,t]*(C1[g]*exp(b[1]+b[2]+b[3]*t))/C3[g]^2/(1-tmp1)+
        n01[g,t]*h2[t]/(1-h2[t]*C3[g])+n00[g,t]/(1-tmp2)*(-exp(b[1]+b[3]*t)*(C2[g]-h2[t]*C1[g])*h2[t]/(1-h2[t]*C3[g])^2)
      }
  }
  return(gr_res)
}


#' Generate simulation data
#' @return \item{n11}{A G by Tn matrix with n11[i,j] being the count of treated subjests with an event within group i at time j. The number of strata is G=5 and the number of time intervals is Tn=20.} 
#' \item{n10}{A G by Tn matrix with n10[i,j] being the count of treated subjests without an event within group i at time j.}
#' \item{n01}{A G by Tn matrix with n01[i,j] being the count of untreated subjests with an event within group i at time j.}
#' \item{n00}{A G by Tn matrix with n00[i,j] being the count of untreated subjests without an event within group i at time j.}
#' @details Besides n11, n10, n01, n00, this function also returns some other simulation paramters, including C1, C2, C3, h2. See Ji et al. (2017) for more details.
#' @references Ji X, Small DS, Leonard CE, Hennessy S (2017). The Trend-in-trend Research Design for Causal Inference. Epidemiology. 28(4), 529–536. 
#' @export
GenData<-function(){
  G<-5
  Tn<-20
  N<-100000
  alpha_0<- -13  #parameters that determine h1 and h2
  #alpha_1<- c(1,0.5,0.1,0.1,0.1)
  alpha_1<- c(1,0.5,0.5,0.5,0.5) # stronger trend
  alpha_2<- 0.9
  alpha_3<-0.008
  x1<-rnorm(N,mean =2, sd = 1)
  x2<-rnorm(N, mean = 2, sd= 1)
  x3<-rbinom(N, 1, 0.8)
  x4<-rbinom(N, 1, 0.2)
  x5<-rbinom(N, 1, 0.1)
  x<-cbind(x1,x2,x3,x4,x5)
  h2<-rep(NA,Tn)
  for(t in 1:15)
    h2[t]<-exp(alpha_2*t-0.5/16*t^2)/10000
  for(t in 16:20)
    h2[t]<-exp(1.8+alpha_2*t-0.6/15.7*t^2)/10000
  h1<-exp(alpha_0 + x%*%alpha_1)*10000
  h1[h1>1/max(h2)]<-1/max(h2)
  x<-x[order(h1),]
  h1<-sort(h1)

  #simulate n0,n1
  Z<-matrix(0,N,Tn)
  for(t in 1:Tn){
    Z[,t]<-rbinom(N,1,h1*h2[t])
  }
  # parameters that determine n
  beta_0<- -4
  beta_1<- log(2.5)
  z<-matrix(0,N,Tn)
  beta_2<- 0.001
  beta_3<-0.1*c(1,0.5,0.1,0.1,0.1)
  #simulate C1 C2 C3 for each group
  C1<-rep(NA,G)
  C2<-rep(NA,G)
  C3<-rep(NA,G)
  for(i in 1:G){
    C1[i]<-mean(exp(x[((i-1)*(N/G)+1):(i*(N/G)),]%*%beta_3)*h1[((i-1)*(N/G)+1):(i*(N/G))])
    C2[i]<-mean(exp(x[((i-1)*(N/G)+1):(i*(N/G)),] %*%beta_3))
    C3[i]<-mean(h1[((i-1)*(N/G)+1):(i*(N/G))])
  }
  #calculate n1, n0 for each group
  n1<-matrix(NA,nrow = G, ncol=Tn)
  n0<-matrix(NA,nrow = G, ncol=Tn)
  for(i in 1:G)
    for(t in 1:Tn){
      n1[i,t]<-sum(Z[((i-1)*(N/G)+1):(i*(N/G)),t])
      n0[i,t]<-N/G-n1[i,t]
    }


  #simulate n11,n10,n01,n00
  n11<-matrix(NA, nrow = G, ncol=Tn)
  n10<-matrix(NA, nrow = G, ncol=Tn)
  n01<-matrix(NA, nrow = G, ncol=Tn)
  n00<-matrix(NA, nrow = G, ncol=Tn)
  for(i in 1:G)
    for(j in 1:Tn){
      n11[i,j]<-rbinom(1,size=n1[i,j],prob=exp(beta_0+beta_1+beta_2*j)*C1[i]/C3[i])
      if(is.na(n11[i,j]))
        print(c(exp(beta_0+beta_1+beta_2*j)*C1[i]/C3[i], n1[i,j]))
      n10[i,j]<-n1[i,j]-n11[i,j]
      n01[i,j]<-rbinom(1,size=n0[i,j],prob=exp(beta_0+beta_2*j)*(C2[i]-C1[i]*h2[j])/(1-C3[i]*h2[j]))
      n00[i,j]<-n0[i,j]-n01[i,j]
    }
  return(list(n11=n11,n10=n10,n01=n01,n00=n00,C1=C1,C2=C2,C3=C3,h2=h2))
}



#######POWER FUNCTION################################

#' Power calculation in trend-in-trend design 
#' @description Monte Carlo power calculation for trend-in-trend design.
#' @param N Sample Size.
#' @param time Number of time points.
#' @param G Number of CPE strata.
#' @param cstat Value of the c-statistic.
#' @param alpha_t A scaler that qunatifies the trend in exposure prevalence.
#' @param beta_0 Intercept of the outcome model.
#' @param h1.OR A given odds ratio.
#' @param nrep Number of Monte Carlo replicates.
#' @return \item{power}{Power of detecting the given Odds Ratio.}
#' @importFrom stats glm predict rbinom
#' @import rms
#' @import pROC
#' @examples
#' \donttest{
#' set.seed(123)
#' ttpower(N=10000,time=10,G=10,cstat=0.75,alpha_t= 0.4,beta_0=-4.3,h1.OR=1.5,nrep=50)
#' }
#' @references Ertefaie A, Small DS, Ji X, Leonard C, Hennessy S (2018). Statistical Power for Trend-in-trend Design. Epidemiology. 29(3), e21–e23. 
#' @export
ttpower<-function(N,time,G,cstat,alpha_t,beta_0,h1.OR,nrep){


expit<-function(x){
	exp(x)/(1+exp(x))
	}

N<-N
Tn<-time
G<-G

alpha_0<- -8
################################################################
################################################################
################################################################
cstat <-cstat

bet.vec<-seq(0,4,0.1)
alpha_2 <-alpha_t
alpha_0<- -8
c<-NULL
y<-matrix(0,N,Tn)
p<-matrix(0,N,Tn)
for( ii in 1:length(bet.vec)){
	bet<-bet.vec[ii]
x<-rnorm(N)

for(t in 1:Tn){
 p[,t]<-exp(alpha_0 + x*bet + alpha_2*t)
 for(i in 1:N)
  if(p[i,t]>1)
   p[i,t]=1
 y[,t]<- rbinom(N,1,p[,t])
}

d <- as.numeric(rowSums(y)>0) 
mod<-glm(d ~x,family="binomial")
predpr<-predict(mod,type="response")

c[ii]<-auc(roc(d ~predpr,plot=FALSE))

}
cbind(bet.vec,c)

alpha_x <-mean(bet.vec[round(c,2)%in%c(cstat-0.01, cstat, cstat +0.01)])




################################################################
################################################################
################################################################




n11<<-NULL
n10<<-NULL
n01<<-NULL
n00<<-NULL
#set.seed(302)
###alpha_0<- -10.5
alpha_0<- alpha_0
alpha_1<- alpha_x
alpha_2<- alpha_t
x<-rnorm(N,mean =-0.0, sd = 1)
u<-rnorm(N,mean =0, sd = 1)

y<-matrix(0,N,Tn)
p<-matrix(0,N,Tn)


for(t in 1:Tn){
 p[,t]<-exp(alpha_0 + x*alpha_1 +alpha_2*t)
 for(i in 1:N)
  if(p[i,t]>1)
   p[i,t]=1
 y[,t]<- rbinom(N,1,p[,t])
}


expo_prev <- colSums(y)/N
expo_prev 


GenOutcome<-function(beta_1=log(1.25)){
 beta_0<- beta_0#-4.3
 beta_1<- beta_1
 z<-matrix(0,N,Tn)
 beta_2<- 0.00
 beta_3<-0
 psum<-rep(NA,Tn)
 zsum<-rep(NA,Tn)
 for(t in 1:Tn){
  p[,t]<-exp(beta_0+beta_1*y[,t]+beta_2*t/3+x*beta_3)#/(1+exp(beta_0+beta_1*y[,t]+beta_2*t+x*beta_3+ u*1.0))
  psum[t]<-sum(p[,t])
  z[,t]<- rbinom(N,1,p[,t]*1) ##########
  zsum[t]<-sum(z[,t])
 }
 h1<-exp(alpha_0 + x*alpha_1 )
 gammax<-exp(x*beta_3)

 outcome_rate<-colSums(z)/N
 outcome_rate 
#plot(outcome_rate,type="o", pch=21,col="blue", xlab="Calendar Quarter", ylab = "Outcome Frequency")

 ps_y <- as.numeric(rowSums(y)>0) # instead 0, use a larger integer be the threshold for ps_y=1 ??
 ps_model_2<-lrm(ps_y~x)
# ps_2<-expit(ps_model_2$coeff[1]+ps_model_2$coeff[2]*x)


 ps_model_2<-lrm(ps_y~1)
 ps_2<-expit(ps_model_2$coeff[1]+x*0)

#auc(roc(ps_y ~ ps_2))


 data_2<-data.frame(y,z,ps_2,h1,gammax)
 ordered_data_2<-data_2[with(data_2, order(ps_2)),]

 C1_2<-rep(NA,G)
 C2_2<-rep(NA,G)
 C3_2<-rep(NA,G)
 for(i in 1:G){
  C1_2[i]<-mean(ordered_data_2$gammax[((i-1)*(N/G)+1):(i*(N/G))]*ordered_data_2$h1[((i-1)*(N/G)+1):(i*(N/G))])
  C2_2[i]<-mean(ordered_data_2$gammax[((i-1)*(N/G)+1):(i*(N/G))])
  C3_2[i]<-mean(ordered_data_2$h1[((i-1)*(N/G)+1):(i*(N/G))])
 }


 #expo prev/ outc freq for every subgroup per calendar quarter

 n11_2<-matrix(NA, nrow = G, ncol=Tn)
 n10_2<-matrix(NA, nrow = G, ncol=Tn)
 n01_2<-matrix(NA, nrow = G, ncol=Tn)
 n00_2<-matrix(NA, nrow = G, ncol=Tn)

 size<-N/G
 for(i in 1:G){
  for(j in 1:Tn){
    temp_y2<-ordered_data_2[c(((i-1)*size+1) : (i*size)) ,j] # treatment
    temp_z2<-ordered_data_2[c(((i-1)*size+1) : (i*size)) ,j+Tn] #outcome
   
    n11_2[i,j] <-sum(temp_y2[temp_z2>0])
    n10_2[i,j] <-sum(temp_y2)-n11_2[i,j]
    n01_2[i,j] <-sum(temp_z2) - n11_2[i,j]
    n00_2[i,j] <- size -sum(temp_y2)-n01_2[i,j]
  } 
 }
 names<- rep(NA,G)
 for(i in 1:G){names[i]<-paste("Quintile",i)}
 plotcol<-rainbow(G)
# par(mar=c(4, 4, 2, 2))
# par(oma=c(1,1,1,1))

# plot((n11_2[G,]+n01_2[G,])/size,type="l",lwd =2, lty = G, main = "Simulated trends in the stratified population " ,xlab="Calendar Quarter", ylab = "Exposure Prevalence",ylim=c(0.015,0.04),cex.lab=1.2)

# for (i in 1:(G-1)){
#  lines((n11_2[i,]+n01_2[i,])/size,type = "l",lwd =2, lty = i, pch=21)
# }
 ##Estimatation of h2(t)
 h2<-rep(NA,Tn)
 h2[1]<-0.001
 for(i in 2:Tn){
  h2[i]<-h2[i-1]*sum(n11_2[,i]+n10_2[,i])/sum(n11_2[,i]+n10_2[,i]+n01_2[,i]+n00_2[,i])*sum(n11_2[,i-1]+n10_2[,i-1]+n01_2[,i-1]+n00_2[,i-1])/sum(n11_2[,i-1]+n10_2[,i-1])
 }
 return(list(n11_2,n10_2,n01_2,n00_2,h2,C1_2,C2_2,C3_2))
}


LL_p<-function(para,para.const){
  b<-rep(NA,3)
  b[1]<-para[1]
  b[2]<-para[2]
  b[3]<-para[3]
#  C1<-para[4:8]
#  C2<-para[9:13]
#  C3<-para[14:18]
  n11<-para.const[[1]]
  n10<-para.const[[2]]
  n01<-para.const[[3]]
  n00<-para.const[[4]]
  h2<-para.const[[5]]
  G<-para.const[[6]]
  Tn<-para.const[[7]]
  loglik<-0
  for(g in 1:G)
    for(t in 1:Tn){
      if(C1[g]<=0)
      {
        C1[g]<-exp(-10)
      }
      if(C3[g]<=0)
      {
        C3[g]<-exp(-10)
      }
      if((C2[g] - C1[g] * h2[t])<0)
      {
        C2[g]<-C1[g] * h2[t]+exp(-10)
      }
      if((1 - C3[g] * h2[t])<0)
      {
        C3[g]<-1/h2[t]-exp(-10)
      }
      if(exp(b[1]+b[3]*t)*(C2[g]-C1[g]*h2[t])/(1-C3[g]*h2[t])>1)
        b[1]<-log(1/(exp(b[3]*t)*(C2[g]-C1[g]*h2[t])/(1-C3[g]*h2[t])))-exp(-10)
      if(exp(b[1]+b[2]+b[3]*t)*C1[g]/C3[g]>1)
        b[1]<-log(1/(exp(b[2]+b[3]*t)*C1[g]/C3[g]))-exp(-10)

      loglik<-loglik+n11[g,t]*(b[1]+b[2]+b[3]*t) + n11[g,t]*(log(C1[g])-log(C3[g]))+
        n10[g,t]*log(1-exp(b[1]+b[2]+b[3]*t)*C1[g]/C3[g])+
        n01[g,t]*(b[1]+b[3]*t)+n01[g,t]*(log(C2[g]-C1[g]*h2[t])-log(1-C3[g]*h2[t]))+
        n00[g,t]*log(1-exp(b[1]+b[3]*t)*(C2[g]-C1[g]*h2[t])/(1-C3[g]*h2[t]))
      
    }
  loglik
}


OR_p<-function(n11,n10,n01,n00){
  G<-dim(n11)[1]
  Tn<-dim(n11)[2]
  n11<-n11
  n10<-n10
  n01<-n01
  n00<-n00
  for(i in 1:G)
    for(j in 1:Tn)
      if(n11[i,j]+n10[i,j]==0)
        n10[i,j]<-1
  h2<-rep(NA,Tn)
  h2[1]<-0.001
  for(i in 2:Tn)
    h2[i]<-h2[i-1]*sum(n11[,i]+n10[,i])/sum(n11[,i]+n10[,i]+n01[,i]+n00[,i])*sum(n11[,i-1]+n10[,i-1]+n01[,i-1]+n00[,i-1])/sum(n11[,i-1]+n10[,i-1])

  b0<-c(-4.3, beta_1,0.00)
  start<-b0
  mle = optim(start,fn=LL_p,gr=NULL,para.const=list(n11,n10,n01,n00,h2,G,Tn,C1,C2,C3),control=list("fnscale"=-1, maxit=10000), method = "Nelder-Mead",hessian=TRUE)

  b1<-mle$par[2]
  sd<-sqrt(abs(diag(solve(-mle$hessian))[2]))
  converge<-mle$converge
  if(converge!=0)
    output<-"The program fails to converge"
  if(converge ==0){
  }
  list(exp(b1),exp(b1-1.96*sd*sqrt(1)),exp(b1+1.96*sd*sqrt(1)))
}



count1<-count2<-ignore<-0
nrep<-nrep
est<-NULL
for(irep in 1:nrep){
      beta_1=log(h1.OR)
    data<-GenOutcome(beta_1)
     n11 <- data[[1]]
     n10 <- data[[2]]
     n01 <- data[[3]]
     n00 <- data[[4]]
     h2 <- data[[5]]
     C1 <- data[[6]]
     C2 <- data[[7]]
     C3 <-  data[[8]]


try(results2 <- OR_p(n11,n10,n01,n00) )
results2
if(results2[[2]]>1| results2[[3]]<1) {count2<-count2+1} 
#if(results2[[3]]>10) ignore<-ignore+1
est[irep]<-results2[[1]]


if(irep%%200==0) print(irep)

}

list(power=count2/(nrep-ignore))

}


#' Finding a detectable odds Ratio with a given power 
#' @description Monte Carlo power calculation for a trend-in-trend design.
#' @param N Sample Size.
#' @param time Number of time points.
#' @param G Number of CPE strata.
#' @param cstat Value of the c-statistic.
#' @param alpha_t A scaler that qunatifies the trend in exposure prevalence.
#' @param beta_0 Intercept of the outcome model.
#' @param power A given power.
#' @param nrep Number of Monte Carlo replicates.
#' @param OR.vec A vector of odds Ratios.
#' @return \item{Power}{A vector of calculated powers for a given OR.vec}
#' \item{OR.vec}{A vector of odds Ratios}
#' \item{DetectDifference }{A detectable difference for a given power value}
#' @references Ertefaie, A., Small, D., Ji, X., Leonard, C., Hennessy, S. (2018). Statistical Power for Trend-in-trend Design. Epidemiology 29(3), e21.
#' @examples 
#' \donttest{
#'set.seed(123)
#'ttdetect(N=10000,time=10,G=10,cstat=0.75,alpha_t= 0.4,beta_0=-4.3,
#'         power=0.80,nrep=50, OR.vec=c(1.9,2.0,2.1,2.2))
#'} 

#' @export
ttdetect<-function(N,time,G,cstat,alpha_t,beta_0,power,nrep, OR.vec){
	pw<-NULL
	detectd<-"NA"
	for(iii in 1:length(OR.vec)){
	temp<-	OR.vec[iii]
	pw[iii]<-ttpower(N=N,time=time,G=G,cstat= cstat,alpha_t= alpha_t,beta_0=beta_0,h1.OR=temp,nrep=nrep)$power
	}

    if(sum(power<=round(pw,2))>0) detectd<-OR.vec[power<=round(pw,2)][1] else detectd<-"The specified power requires larger OR"
	list(Power=pw, OR.vec= OR.vec, DetectDifference=detectd)

}



