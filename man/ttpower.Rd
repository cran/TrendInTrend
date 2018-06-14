 % File src/library/base/man/load.Rd
\name{ttpower}
\alias{ttpower} 
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Power calculation in trend-in-trend design 
}
\description{
Monte Carlo power calculation for a trend-in-trend design


}
\usage{
ttpower(N,time,G,cstat,alpha_t,beta_0,h1.OR,nrep)
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
  \item{h1.OR}{
A given odds ratio.
}
  \item{nrep}{
Number of Monte Carlo replicates.
}
}
\value{
%%  ~Describe the value returned
Power of detecting the given Odds Ratio
%% ...
}
\references{
Ertefaie, A., Small, D. S., Ji, X., Leonard, C., & Hennessy, S. (2018) Statistical Power for Trend-in-trend Design. Epidemiology (Cambridge, Mass.), 29(3), 21-23.
}
\examples{
\dontrun{
set.seed(123)
ttpower(N=10000,time=10,G=10,cstat=0.75,alpha_t= 0.4,beta_0=-4.3,h1.OR=1.5,nrep=50)
  }
   } 

