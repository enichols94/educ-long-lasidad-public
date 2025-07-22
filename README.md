# Code to accompany: "The association between education and cognitive decline: An exploration of unexpected findings in a nationally representative longitudinal study in India"

Code includes data processing and all analytic code corresponding to analysis for the above paper. Files are organized approximately as such: 

- 010_dataprep.R: data preparation
- 011_mort_weights.R: estimate IPW weights for attrition
- 012_pe_weights.R: estimate IPW weights to make the refresher and returner samples comparable to calculate practice effects
- 020_descriptivetable.R: create Table 1
- 021_plottrajectories.R: plot raw trajectories for Figure 1
- 030_longmodels.R: estimate models included across analysis - adjusted models, stratified models, models for different cognitive outcomes, models accounting for mortality
- 040_regressiontomean.R: descriptive plot for regression to the mean, incorporate simulation results (saved in 061_simulation_rtm.R)
- 050_practiceeffects.R: estimate practice effects
- 060_simulation_survival.R: simualtion to evaluate bias due to selective survival
- 061_simulation_rtm.R: simulation to evaluate scenarios with different correlation between wave 1 and change

  
