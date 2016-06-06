#' An Odds Ratio Estimation Function
#' @description estimate causal odds ratio given trends in exposure prevalence and outcome frequencies of stratified data.
#' @param n11 A G by T matrix with n11[i,j] being the counts of positive outcomes among treated subjests within group i at time j;
#' @param n10 A G by T matrix with n10[i,j] being the counts of negative outcomes among treated subjests within group i at time j;
#' @param n01 A G by T matrix with n01[i,j] being the counts of positive outcomes among control subjests within group i at time j;
#' @param n00 A G by T matrix with n00[i,j] being the counts of negative outcomes among control subjests within group i at time j.
#' @return ORs ratio and 95\% conference interval
#' @importFrom stats optim rnorm rbinom
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
OR<-function(n11,n10,n01,n00){
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

  b0<-c(-4,0,0)
  G0<-rep(c(0.5,1,0.2),rep(5,3))
  start<-c(b0,G0)
  mle = optim(start,fn=LL,gr=NULL,para.const=list(n11,n10,n01,n00,h2,G,Tn),control=list("fnscale"=-1, maxit=10000), method = "Nelder-Mead",hessian=TRUE)
  b1<-mle$par[2]
  sd<-sqrt(abs(diag(solve(-mle$hessian))[2]))
  converge<-mle$converge
  if(converge!=0)
    output<-"The program fails to converge"
  if(converge ==0){
    output<-list(paste("Odds Ratio:",exp(b1)), paste("Confidence Interval:",exp(b1-1.96*sd), exp(b1+1.96*sd)))
  }
  output
}

LL<-function(para,para.const){
  b<-rep(NA,3)
  b[1]<-para[1]
  b[2]<-para[2]
  b[3]<-para[3]
  C1<-para[4:8]
  C2<-para[9:13]
  C3<-para[14:18]
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

#' Generate simulation data
#' @export
GenData<-function(){
  G<-5
  T<-20
  N<-100000
  alpha_0<- -13  #parameters that determine h1 and h2
  alpha_1<- c(1,0.5,0.1,0.1,0.1)
  alpha_2<- 0.9
  alpha_3<-0.008
  x1<-rnorm(N,mean =2, sd = 1)
  x2<-rnorm(N, mean = 2, sd= 1)
  x3<-rbinom(N, 1, 0.8)
  x4<-rbinom(N, 1, 0.2)
  x5<-rbinom(N, 1, 0.1)
  x<-cbind(x1,x2,x3,x4,x5)
  h2<-rep(NA,T)
  for(t in 1:15)
    h2[t]<-exp(alpha_2*t-0.5/16*t^2)/10000
  for(t in 16:20)
    h2[t]<-exp(1.8+alpha_2*t-0.6/15.7*t^2)/10000
  h1<-exp(alpha_0 + x%*%alpha_1)*10000
  h1[h1>1/max(h2)]<-1/max(h2)
  x<-x[order(h1),]
  h1<-sort(h1)

  #simulate n0,n1
  Z<-matrix(0,N,T)
  for(t in 1:T){
    Z[,t]<-rbinom(N,1,h1*h2[t])
  }
  # parameters that determine n
  beta_0<- -4
  beta_1<- 0#log(2.5)
  z<-matrix(0,N,T)
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
  n1<-matrix(NA,nrow = G, ncol=T)
  n0<-matrix(NA,nrow = G, ncol=T)
  for(i in 1:G)
    for(t in 1:T){
      n1[i,t]<-sum(Z[((i-1)*(N/G)+1):(i*(N/G)),t])
      n0[i,t]<-N/G-n1[i,t]
    }


  #simulate n11,n10,n01,n00
  n11<-matrix(NA, nrow = G, ncol=T)
  n10<-matrix(NA, nrow = G, ncol=T)
  n01<-matrix(NA, nrow = G, ncol=T)
  n00<-matrix(NA, nrow = G, ncol=T)
  for(i in 1:G)
    for(j in 1:T){
      n11[i,j]<-rbinom(1,size=n1[i,j],prob=exp(beta_0+beta_1+beta_2*j)*C1[i]/C3[i])
      if(is.na(n11[i,j]))
        print(c(exp(beta_0+beta_1+beta_2*j)*C1[i]/C3[i], n1[i,j]))
      n10[i,j]<-n1[i,j]-n11[i,j]
      n01[i,j]<-rbinom(1,size=n0[i,j],prob=exp(beta_0+beta_2*j)*(C2[i]-C1[i]*h2[j])/(1-C3[i]*h2[j]))
      n00[i,j]<-n0[i,j]-n01[i,j]
    }
  return(list(n11,n10,n01,n00,C1,C2,C3,h2))
}

