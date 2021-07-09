#' Pipe operator
#' 
#' Define magrittr's Pipe operator "%>%" locally for simplified programming.
#'
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' 
#' @return
#' Function call rhs(lhs).
#' 
`%>%` <- magrittr::`%>%`