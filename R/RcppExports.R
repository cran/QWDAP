# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

stepOne <- function(findIn, p, n, sigma, tolerance, Ftrace, criteria, Y, X1, X0, k, SST) {
    .Call('_QWDAP_stepOne', PACKAGE = 'QWDAP', findIn, p, n, sigma, tolerance, Ftrace, criteria, Y, X1, X0, k, SST)
}

qwalkRcpp <- function(edges, startindex, lens, scals, flag, getfloat, multiple) {
    .Call('_QWDAP_qwalkRcpp', PACKAGE = 'QWDAP', edges, startindex, lens, scals, flag, getfloat, multiple)
}

