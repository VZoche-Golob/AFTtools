## Additional tools for Parametric Survival Models (a.k.a. Accelerated Failure Time Models)

The [R](https://www.r-project.org/) package **AFTtools** provides additional tools
to work with parametric survival models fitted by [`survival::survreg`](https://cran.r-project.org/web/packages/survival/).  

To select the most appropriate distribution for the model, several models with 
different distributions can be fit with only one call of `compare_PSMdist()`. Parallel 
computation can be used to fasten the model estimations. The result can be printed 
as a nice table by `show_comparison()`.  

Population average predictions (survival and quantiles) from parametric survival 
models with frailties can be made using `predict_survreg()`. Bootstrap confidence 
intervals of the predictions are also provided.  

The functions `AFTplot()`, `CSRplot()`, `Iplot()`, and `NPplot()` produce diagnostic 
plots.  

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.46383.svg)](http://dx.doi.org/10.5281/zenodo.46383) (Release version 0.2.1, see [NEWS](NEWS.md) for version changes.)   

***

#### Note:  

`NPplot()` depends on [**interval**](https://cran.r-project.org/web/packages/interval/index.html) which requires the package [**Icens**](http://www.bioconductor.org/packages/release/bioc/html/Icens.html) from [Bioconductor](http://www.bioconductor.org/).