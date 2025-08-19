library(data.table)




get_enrol_time <- function(N = 2500, lambda = 1.52,
                               rho = function(t) pmin(t/360, 1)){
  
  c(0, poisson::nhpp.event.times(lambda, N - 1, rho))
}



