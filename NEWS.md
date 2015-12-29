## Version log of **AFTtools**



#### Version 0.2.1
finished: 2015-12-29  

##### Fixes
* `predict_survreg` : Fixed the transfer of `strata = NULL` to `boot::boot`.  
* `AFTplot` : Title on all plot pages.  
* `NPplot` : Fixed handling of interaction terms in model formulas.  

##### Improvements
* `AFTplot` : Added quantile regression line to facilitate interpretation.  



### Version 0.2
finished: 2015-11-30

##### Fixes
* `CSRplot` : Calculation of the Cox-Snell-Residuals corrected according to Farrington (2000)  
* `create_int2Surv` : Right open intervals corrected  
* `predict_survreg` : Fixed bug that resulted in confidence intervals not covering the point estimates for distributions other than Weibull   