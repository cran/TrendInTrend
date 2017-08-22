 % File src/library/base/man/load.Rd
\name{ttdetect}
\alias{ttdetect} 
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Finding a detectable odds Ratio with a given power 
}
\description{
Monte Carlo power calculation for a trend-in-trend design


}
\usage{
ttdetect(N,time,G,cstat,alpha_t,beta_0,power,nrep, OR.vec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
Sample size.
}
  \item{time}{
Number of time points.
}
  \item{G}{
Number of CPE strata
}
  \item{cstat}{
Value of the c-statistic.
}
  \item{alpha_t}{
A scaler that qunatifies the trend in exposure prevalence.
}
  \item{beta_0}{
Intercept of the outcome model.
}
  \item{power}{
A given power.
}
  \item{nrep}{
Number of Monte Carlo replicates.
}
  \item{OR.vec}{
A vector of odds Ratios
}

}
\value{
%%  ~Describe the value returned
An object 'ttdetect' is a list of containing the following components:
  \item{Power }{A vector of calculated powers for a given OR.vec }
  \item{OR.vec }{A vector of odds Ratios}
  \item{DetectDifference }{A detectable difference for a given power value}
%% ...
}


\examples{
\dontrun{
ttdetect(N=10000,time=10,G=10,cstat=0.75,alpha_t= 0.4,beta_0=-4.3,
power=0.80,nrep=50, OR.vec=c(1.9,2.0,2.1,2.2))
  } 
   } 

