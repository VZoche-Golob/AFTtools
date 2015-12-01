## Version log of **AFTtools**



#### Version 0.2.1 (under development)


##### Fixes
* `predict_survreg` : Fixed the transfer of `strata = NULL` to `boot::boot`.  



### Version 0.2
2015-11-30

##### Fixes
* `CSRplot` : Calculation of the Cox-Snell-Residuals corrected according to Farrington (2000)  
* `create_int2Surv` : Right open intervals corrected  
* `predict_survreg` : Fixed bug that resulted in confidence intervals not covering the point estimates for distributions other than Weibull  

##### Notes
* Parametric bootstrap samples don’t seem to be appropriate for `MIC`.