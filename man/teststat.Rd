% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrappers.R
\name{teststat}
\alias{teststat}
\title{Calculate test statistic}
\usage{
teststat(df, n_E, n_C, delta, method, better)
}
\arguments{
\item{df}{data frame with variables \mjseqn{x_E} and \mjseqn{x_C}}

\item{n_E}{Sample size in experimental group.}

\item{n_C}{Sample size in control group.}

\item{delta}{Non-inferiority margin.}

\item{method}{Specifies the effect measure/test statistic. One of "RD", "RR", or "OR".}

\item{better}{"high" if higher values of x_E favour the alternative
hypothesis and "low" vice versa.}
}
\value{
A list with the two elements \code{p_max} and \code{p_vec}.
\code{p_max} is the maximum p-value and most likely serves as "the one" p-value.
\code{p_vec} is a named vector. The names indicate the true proportion pairs
\mjseqn{(p_E, p_C)} with \mjseqn{e(p_E, p_C) = \delta} that underly the calculation of
the p-values. It can be used for plotting the p-value versus the true proportions.
}
\description{
\loadmathjax
\code{teststat} takes a data frame with variables \code{x_E} and \code{x_C}
and adds the variable \code{stat} with the value of the test statistic
specified by \code{n_E}, \code{n_C}, \code{method}, \code{delta} and \code{better}.
}
\details{
If higher values of \mjseqn{x_E} favour the alternative hypothesis (\code{better = "high"}), we are interested
in testing the null hypothesis
\mjsdeqn{H_0: e(p_E, p_C) \le \delta ,}
where \mjseqn{e} is one of the effect measures risk difference (\code{method = "RD"}),
risk ratio (\code{method = "RR"}), or odds ratio (\code{method = "OR"}).
The test statistic for risk difference is
\mjsdeqn{T_{\mbox{RD}, \delta}(x_E, x_C) = \frac{\hat{p_E} - \hat p_C - \delta}{\sqrt{\frac{\tilde p_E(1 - \tilde p_E)}{n_E} + \frac{\tilde p_C(1 - \tilde p_C)}{n_C}}},}
where \mjseqn{\tilde p_C = \tilde p_C(x_E, x_C)} is the MLE of \mjseqn{p_C} and
\mjseqn{\tilde p_E = \tilde p_C + \delta} is the MLE of \mjseqn{p_E} under \mjseqn{p_E - p_C = \delta}.
High values of \mjseqn{T_{\mbox{RD}, \delta}} favour the alternative hypothesis.
The test statistic for risk ratio is
\mjsdeqn{T_{\mbox{RR}, \delta}(x_E, x_C) = \frac{\hat p_E - \delta \cdot \hat p_C}{\sqrt{\frac{\tilde p_E(1 - \tilde p_E)}{n_E} + \delta^2\frac{\tilde p_C(1 - \tilde p_C)}{n_C}}},}
where \mjseqn{\tilde p_C = \tilde p_C(x_E, x_C)} is the MLE of \mjseqn{p_C} and
\mjseqn{\tilde p_E = \tilde p_C + \delta} is the MLE of \mjseqn{p_E} under \mjseqn{p_E / p_C = \delta}.
High values of \mjseqn{T_{\mbox{RR}, \delta}} favour the alternative hypothesis.
The test statistic for Odds Ratio
\mjsdeqn{ T_{\mbox{OR}, \delta} = 1-(1 - F_{\mbox{ncHg}(X_E+X_C, n_E, n_C, \delta)}(x_E-1)) }
is based on Fisher's non-central hypergeometric distribution with density
\mjsdeqn{ f_{\mbox{ncHg}(s, n_E, n_C, \delta)}(k) = \frac{\binom{n_E}{k}\cdot \binom{n_C}{s-k}\cdot \delta^k}{\sum\limits_{l \in A_{s, n_E, n_C}} \binom{n_E}{l}\cdot \binom{n_C}{s-l}\cdot \delta^l}, }
where \mjseqn{A_{s, n_E, n_C} = \{\max(0, s-n_C), \dots, \min(n_E, s)\}}.
The density is zero if \mjseqn{k < \max(0, s-n_C)} or \mjseqn{k > \min(n_E, s)}.
High values of \mjseqn{T_{\mbox{OR}, \delta}} favour the alternative hypothesis (due to "1-...").
}
\examples{
n_E <- 10
n_C <- 11
df <- expand.grid(
  x_E = 0:n_E,
  x_C = 0:n_C
)
teststat(
  df = df,
  n_E = n_E,
  n_C = n_C,
  method = "RD",
  delta = -0.1,
  better = "high"
)
}
