% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boschloo-functions.R
\name{samplesize_Wang}
\alias{samplesize_Wang}
\title{Calculate approximate sample size for approximate test}
\usage{
samplesize_Wang(p_EA, p_CA, gamma, alpha, beta, r)
}
\arguments{
\item{p_EA}{assumed response rate in the experimental group under the alternative}

\item{p_CA}{assumed response rate in the control group under the alternative}

\item{gamma}{TODO}

\item{alpha}{maximal type I error rate}

\item{beta}{maximal type II error rate}

\item{r}{allocation ratio}
}
\value{
Sample sizes per group (n_C, n_E).
}
\description{
alculate approximate sample size for approximate test for specified
level alpha, power, allocation ratio r = n_E/n_C, true rates p_CA, p_EA and
OR-NI.margin gamma.
}
