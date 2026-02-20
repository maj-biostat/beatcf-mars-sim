



#' Inhomogeneous Poisson process (Ogata’s algorithm)
#' 
#' Generate event times up to specified number of events configured to 
#' emulate an enrolment process where there is a linear ramp up to a maximum
#' rate parameter after which the process functions as a homogenous PP with 
#' fixed rate.
#' 
#' IPP has same properties as poisson process, except for the fact that its rate 
#' is a function of time, i.e., \lambda=\lambda(t)
#' 
#' For example, a football team may score goals at a higher rate at the end of 
#' a game than at the beginning if the opposing team tires more quickly. 
#' Similarly, a clinical trial may have slow recruitment at the beginning and
#' gradually ramp up to a stable rate (generally this is a bit artificial but
#' is probably better than just assuming a constant enrolment rate from the 
#' start). 
#' 
#' Define Inhomogeneous Poisson process as a counting process: {N(t), t >= 0}, 
#' so that N has integer values that never decrease over time, but jump up at 
#' random times, and we specify the following 4 conditions:
#' 
#' N(0) = 0
#' increments are independent (Markov) but not stationary
#' P(N(t + h) − N(t) = 1) = λ(t)h + o(h)
#' P(N(t + h) − N(t) > 1) = o(h)
#'
#' The number of arrivals in any interval is a Poisson random variable but 
#' the parameter can depend on the location of the interval.
#' 
#' N(t+s)−N(t) \sim \text{Poisson}(\int_t^{t+s} \lambda(a) da)
#' 
#' Ogata algorithm involves:
#' 
#' Find a constant \lambda_{max} >= \lambda(t)
#' Simulate a HPP with rate \lambda_{max}
#' Keep each event at time t with with probability \lambda(t)/\lambda_{max}
#' 
#' where \lambda_{max} is the rate at which the enrolment stabilises.
simulate_ipp_fixed_n <- function(
    N,
    lambda = function(t, lambda_inf = 1.52, ramp_up = 90) {
      lambda_inf * pmin(t/ramp_up, 1)
    },
    lambda_max,
    ramp_up_period
) {
  t <- 0
  events <- numeric(0)
  while (length(events) < N) {
    # Propose next event from homogeneous PP
    t <- t + rexp(1, rate = lambda_max)
    
    # Accept with probability lambda(t) / lambda_max
    lambda_t <- lambda(t, lambda_max, ramp_up_period)
    if (runif(1) < lambda_t / lambda_max) {
      events <- c(events, t)
    }
  }
  # event times
  events
}



