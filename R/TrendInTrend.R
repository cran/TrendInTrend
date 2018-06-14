#' An Odds Ratio Estimation Function
#' @description estimate causal odds ratio given trends in exposure prevalence and outcome frequencies of stratified data.
#' @param n11 A G by T matrix with n11[i,j] being the counts of positive outcomes among treated subjests within group i at time j;
#' @param n10 A G by T matrix with n10[i,j] being the counts of negative outcomes among treated subjests within group i at time j;
#' @param n01 A G by T matrix with n01[i,j] being the counts of positive outcomes among control subjests within group i at time j;
#' @param n00 A G by T matrix with n00[i,j] being the counts of negative outcomes among control subjests within group i at time j;
#' @param bnull The initial values for the optimization algorithm.
#' @return ORs ratio, 95\% conference interval, log-likelihood value
#' @importFrom stats optim rnorm rbinom 
#' @importFrom stats glm predict
#' @importFrom grDevices rainbow
#' @import pROC
#' @examples
#' \donttest{
#' data <- GenData()
#' n11 <- data[[1]]
#' n10 <- data[[2]]
#' n01 <- data[[3]]
#' n00 <- data[[4]]
#' results <- OR(n11,n10,n01,n00,bnull=c(-4,0,0))
#' }
#' @export
OR<-function(n11,n10,n01,n00,bnull=c(-4,0,0)){
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

  b0<-bnull
  G0<-rep(c(0.5,1,0.2),rep(G,3))
  start<-c(b0,G0)
  mle = optim(start,fn=LL,gr=NULL,para.const=list(n11,n10,n01,n00,h2,G,Tn),control=list("fnscale"=-1, maxit=10000), method = "Nelder-Mead",hessian=TRUE)
  b1<-mle$par[2]
  sd<-sqrt(abs(diag(solve(-mle$hessian))[2]))
  converge<-mle$converge
  LV<- round(LL(mle$par,para.const=list(n11,n10,n01,n00,h2,G,Tn)),2)
  if(converge!=0)
    output<-"The program fails to converge"
  if(converge ==0){
    output<-list(paste("Odds Ratio:",exp(b1)), paste("Confidence Interval:",exp(b1-1.96*sd), exp(b1+1.96*sd)),paste("Log-Likelihood Value:",LV))
  }
  output
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













#######POWER FUNCTION################################

#' A power Function
#' @description estimate power.
#' @return power
#' @importFrom stats glm predict rbinom
#' @import rms
#' @import pROC
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
#par(mar=c(4, 4, 2, 2))
#par(oma=c(1,1,1,1))
#plot(expo_prev,type="l", pch=21,lwd =2, col="black", main="Simulated trend in the entire population",xlab="Calendar Quarter", ylab = "Exposure Prevalence",cex.lab=1.2)



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



#set.seed(123)
#ttpower(N=10000,time=10,G=10,cstat=0.75,alpha_t= 0.4,beta_0=-4.3,h1.OR=1.5,nrep=100)






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

#set.seed(123)
#ttdetect(N=10000,time=10,G=10,cstat=0.75,alpha_t= 0.4,beta_0=-4.3,power=0.80,nrep=10, OR.vec=c(2.0,2.1,2.2))

