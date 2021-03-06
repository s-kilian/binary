% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boschloo-functions.R
\name{samplesize_exact_boschloo}
\alias{samplesize_exact_boschloo}
\title{Calculate exact sample size for Fisher-Boschloo test}
\usage{
samplesize_exact_boschloo(
  p_EA,
  p_CA,
  alpha,
  beta,
  r,
  size_acc = 4,
  alternative = "greater"
)
}
\arguments{
\item{p_EA}{assumed response rate in the experimental group under the alternative}

\item{p_CA}{assumed response rate in the control group under the alternative}

\item{alpha}{maximal type I error rate}

\item{beta}{maximal type II error rate}

\item{r}{allocation ratio}

\item{size_acc}{TODO}

\item{alternative}{character indicating if the alternative is "greater" or "smaller"}
}
\value{
Sample sizes per group (n_C, n_E), nominal alpha and exact power.
}
\description{
Calculate exact sample size for Fisher-Boschloo test and specified
level alpha, power, allocation ratio r = n_E/n_C and true rates p_CA, p_EA.
Accuracy of calculating the critical value can be specified by size_acc.
}
