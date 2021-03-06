% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TrendInTrend.R
\name{OR}
\alias{OR}
\title{An Odds Ratio Estimation Function}
\usage{
OR(n11, n10, n01, n00, bnull = c(-10, 0, 0), n_explore = 10,
  noise_var = c(1, 1, 0.5), n_boot = 50, alpha = 0.05)
}
\arguments{
\item{n11}{A G by Tn matrix with n11[i,j] being the count of treated subjects with an event within group i at time j. The number of strata is G and the number of time intervals is Tn.}

\item{n10}{A G by Tn matrix with n10[i,j] being the count of treated subjects without an event within group i at time j.}

\item{n01}{A G by Tn matrix with n01[i,j] being the count of untreated subjects with an event within group i at time j.}

\item{n00}{A G by Tn matrix with n00[i,j] being the count of untreated subjects without an event within group i at time j.}

\item{bnull}{Initial values for beta0, beta1, beta2 for the optimization algorithm. Default is (-10,0,0). It is suggested the initial value of beta0 be set as a small negative number (-4 or smaller) for the rare outcome model to be computationally stable.}

\item{n_explore}{Number of iterations in the optimization algorithm to stabilize the outputs. Default is 10.}

\item{noise_var}{The optimization algorithm is iterated n_explore times. Results from the previous iteration with added Gaussian noise are set as the starting values for the new iteration. Bigger noise_var indicates larger variance for the Gaussian noise, meaning more exploration during the iterations. Default is (1,1,0.5).}

\item{n_boot}{Number of bootstrap iterations to construct the confidence interval for the estimated odds ratio beta1. Default is 50.}

\item{alpha}{(1-alpha) is the significance level of the confidence interval. Default is 0.05.}
}
\value{
\item{beta}{Maximum likelihood estimators (MLE) for beta0, beta1, beta2. Beta1 is the estimated treatment-event odds ratio. Because we conduct n_explore iterations, the set of parameters that is associated with the highest log likelihood is the output.} 
\item{CI_beta1}{1-alpha confidence interval for beta1.}
\item{ll}{Log likelihood evaluated at the MLE.}
\item{not_identified}{Equals 1 if the MLE is not identifiable or weakly identified. This could happen when there are multiple sets of parameters associated with the highest log likelihood, or the bootstrap confidence interval fails to cover the estimated beta1.}
}
\description{
Estimate causal odds ratio (OR) given trends in exposure prevalence and outcome frequencies of stratified data.
}
\details{
This function estimates the odds ratio parameter beta1 in the subject-specific model in Ji et al. (2017)
\deqn{logit(E[Y(it)|Z(it), G(i), X(it)])=beta0+Z(it)*beta1+t*beta2+X(it)\gamma}
where \eqn{Z(it)} and \eqn{Y(it)} are the binary exposure and outcome variables for individual \eqn{i} at time \eqn{t}. 
There are three caveats regarding the implementation. First, the trend-in-trend design works better when there are substantial exposure trend differences across strata. If the exposure trend is roughly parallel across strata, the method may fail to converge. Second, we recommend running the OR function for multiple starting points to evaluate the stability of the optimization algorithm. Third, the bootstrap confidence interval may have slightly lower coverage probability than the nominal significance level 1-alpha.
}
\examples{
\donttest{
data <- GenData()
n11 <- data[[1]]
n10 <- data[[2]]
n01 <- data[[3]]
n00 <- data[[4]]
results <- OR(n11,n10,n01,n00)
}
}
\references{
Ji X, Small DS, Leonard CE, Hennessy S (2017). The Trend-in-trend Research Design for Causal Inference. Epidemiology 28(4), 529–536.  \cr
Ertefaie, A., Small, D., Ji, X., Leonard, C., Hennessy, S. (2018). Statistical Power for Trend-in-trend Design. Epidemiology 29(3), e21.\cr
Ertefaie, A., Small, D., Leonard, C., Ji, X., Hennessy, S. (2018). Assumptions Underlying the Trend-in-Trend Research Design. Epidemiology  29(6), e52-e53.
}
