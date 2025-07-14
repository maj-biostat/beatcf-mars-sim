# beatcf-mars-sim

Simulation work for beat cf mars.

## Overview

Trial design is a Bayesian, open-label, parallel group, 3-arm, sample-size adaptive, randomised non-inferiority trial nested within BEAT-CF platform.

Three-arm adaptive trial evaluating the relative effect of deferred antibiotic or early stopping relative to the standard of care therapy for pulmonary exacerbations.

+ Standard of care is immediate start of oral antibiotic at clinical discretion with the plan to stop after 14 days, regardless of response.
+ Immediate start of oral antibiotic at clinical discretion with the plan to stop once symptoms are resolved
+ Defer start of oral antibiotic and start based on pre-specified rules

Outcome is binary with failure event corresponding to self report worsening of symptoms at day 4 (relative to day 0), no improvement at day 7 or not recovered at day 14.
Approximate 30% of participants are expected to have failures.

Population is stratified based on lung function in the 12 months prior to the exacerbation and airway colonisation (pseudomonas).
The following distribution has been assumed.

| Group | Lung Function | Pseudomonas Status | Estimated Proportion | 
| ----- | ------------- | ------------------ | -------------------- | 
| A     | ≥ 70%         | Not Colonised      | Approximately 40%    | 
| B     | ≥ 70%         | Colonized          | Approximately 20%    | 
| C     | < 70%         | Not Colonised      | Approximately 15%    | 
| D     | < 70%         | Colonized          | Approximately 25%    | 

The randomisation is 1:1:1, stratified by stratum membership at the time of enrolment.

For the analysis, participants are analysed according to their assigned strategy.
Currently the design considers only the first exacerbation, but a variation may explore a re-randomisation design.
While the primary analysis model will include variables for stratum, age, sex and site, the simulation model only includes stratum membership.

The non-inferiority margin is set to 0.15 in absolute terms on the probability scale.
If a treatment arm shows a 0.99 probability of being non-inferior to the SoC, the strategy is deemed non-inferior.
Inferiority rules are also implemented.

For each of the two test arms:

+ If non-inferiority is demonstrated:
  + Stop enrollment to this arm and consider recommending this as an acceptable alternative to standard treatment
+ If inferiority is demonstrated:
  + Stop enrollment to this arm
  
In the final results, the treatments are ranked based on their respective probability of treatment success.

Analyses start at 300 participants having reached the 14 day follow up for the first exacerbation with subsequent analyses at 375 and 450.



