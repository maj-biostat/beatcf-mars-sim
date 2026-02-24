# Markov models


Inhomogeneous Markov models (NHMM) and semi-Markov models (SMM) both relax the strict assumptions of standard Markov chains by introducing time-dependence, but differ in what that dependency is based on. NHMMs depend on absolute time (calendar time), while SMMs depend on the duration spent in the current state. 
SMMs allow non-exponential, general distribution for holding times, providing greater flexibility in modeling aging or duration-dependent processes. 

+ Time-Varying Probabilities: The transition probability matrix $P(s, t)$ depends on both the starting time and ending time $t$.
+ Non-Stationary: The process does not settle into a constant, long-term behavior.
+ Memoryless Property: Despite changing over time, the future state still depends only on the current state, not the past history.
+ Modeling Approaches:
  + Discrete Time: Uses different transition matrices for different time steps.
  + Continuous Time: Uses time-dependent transition rate matrices $Q(t)$
+ Applications:
  + Resource Availability: Modeling daily or weekly fluctuations in system loads.
  + Reliability Engineering: Modeling non-linear, time-dependent degradation.
  + Population Genetics: Describing genetic drift in finite, changing populations.
  + Vehicle Driving Patterns: Capturing diurnal variations in usage. 


Differences from Homogeneous Models:
While homogeneous Markov chains are simpler to analyze and have constant transition probabilities, inhomogeneous chains allow for more realistic modeling of dynamic systems. 
However, inhomogeneous models are generally more computationally intensive to simulate and estimate, often requiring specialized techniques like adapted Baum-Welch algorithms for hidden Markov scenarios. 

The time-inhomogeneous and semi-Markov case might look similar, but they are different.

There are two 'clocks' in play. 
There is the 'time' clock, which ticks continuously throughout time, and then there is the 'duration' clock, which is reset at every jump. 
If you imagine that $X$ are health events of an individual, what is going one can be explained as follows. 

If $X$ is time-inhomogeneous Markov, the time until the next health event (healthy, ill, dead) is allowed to depend on:

+ the current health state of the individual (health, ill, dead)
+ his age (time). 

If $X$ is (time-homogeneous) semi-Markov, the time until the next health event is allowed to depend on:

+ the current health state of the individual
+ the time since he entered this state of health (duration). 

These two characteristics are of course very different.

Both aspects are combined in the **time-inhomogeneous semi-Markov** case with transition rates $\mu_{jk}(t,u)$ (that depend on time as well as duration). 
Here the time until the next health event is allowed to depend on:

+ the current health state of the individual, 
+ the time since he entered this state of health, and 
+ his age.


## Reminders for censoring and truncation

https://www.geeksforgeeks.org/data-analysis/censoring-and-truncation/
https://blog.stata.com/2016/12/13/understanding-truncation-and-censoring/


Our data are left-truncated when individuals below a threshold are not present in the sample. 
For example, if we want to study the size of certain fish based on the specimens captured with a net, fish smaller than the net grid won’t be present in our sample.

Our data are left-censored at 𝜅 if every individual with a value below 𝜅 is present in the sample, but the actual value is unknown. 
This happens, for example, when we have a measuring instrument that cannot detect values below a certain level.

Left censoring occurs if a participant is entered into the study when the milestone of interest occurred prior to study entry but the age at that milestone is unknown. 
Left truncation occurs when individuals who have already passed the milestone at the time of study recruitment are not included in the study.


Left-censoring occurs
when we cannot observe the time when the event occurred. For obvious reasons if the event is death,
the data can’t be left-censored. A good example is discussed in an ASA paper on survival analysis, “e.g.
[a] study of age at which African children learn a task. Some already knew (left-censored), some learned
during a study (exact), some had not yet learned by end of study (right-censored).”

Truncation is due to sampling bias that only those individuals whose lifetimes lie within a certain interval can be observed.
Left truncation occurs in data collection when observations below a certain threshold are completely excluded from the sample, often because they failed to survive or exist until a specific, later observation point. 
For example, a study on cancer mortality begins in 2025, only including patients alive at that time. 
Patients who were diagnosed and died before 2025 are not included, resulting in a sample biased toward longer survival times.



