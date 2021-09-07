# parTreat

Efficient estimation of treatment effects based on parametric models of the treatment effect, as proposed by Athey et al. (2021).  

To install this package in R, run the following commands:  

```R
library(devtools) 
install_github("michaelpollmann/parTreat")
```



## Usage

If `X` are the outcomes of the control observations, and `Y` are the outcomes of the treated observations, call

```r
eif_additive(X,Y)
```
or
```r
waq(X,Y)
```
for the influence function-based or the weighted average of quantile differences estimators, to efficiently estimate the treatment effect assuming additivity. The partially adaptive asymmetric trimmed mean estimator is available as
```r
atm_diff(X,Y)
```


For a multiplicative model, use
```r
est_eif_log <- eif_additive(log(X),log(Y))
log_to_level(est_eif_log, X, Y)
```
or
```r
est_waq_log <- waq(log(X),log(Y))
log_to_level(est_waq_log, X, Y)
```



## Brief Description

Efficiently estimate treatment effects imposing a parametric assumption on treatment effects of the form Y(1) = h(Y(0), a) where a is the parameter to be estimated.
Currently, three estimators for the additive model Y(1) = Y(0) + tau are implemented, where tau is the treatment effect.
There is also a helper function to estimate an additive specification in logs and translate the estimated effect and standard error to level-effects, such that the functions can be used for the multaplicative model Y(1) = a * Y(0) as well.

The three estimators are:  

- `eif_additive`: influence-function based estimator, which uses the first derivative of the log density
- `waq`: weighted average of quantile differences, which uses weights proportional to (minus) the second derivative of the log density to efficiently weight the differences in quantiles (order statistics) of treated and control
- `atm_diff`: difference in asymmetric trimmed means of `Y` and `X`; where `Y` and `X` are trimmed in the same way, but the trimming percentages for left and right tails can be different, which can be convenient for non-negative data where primarily the right tail is of concern (set argument `max_alpha=0` to not trim the left tail). For computational feasibility with large samples, by default the function searches over the optimal trimming parameters with a `sqrt(n)` grid; for an exhaustive search, set argument `gridsize='n'`.

Under correct specification (of the additive or multiplicative model), the eif and waq estimators are asymptotically equivalent.
They are adaptive, in the sense that if the additive model is correctly specified,   they are just as efficient as parametric estimators that use knowledge of the true distribution of the outcomes (up to shift).



## Note

The `eif_additive` and `waq` functions estimate the density based on the control observations.
This is convenient if the control group is larger than (or at least as large as) the treatment group.
If the treated group is larger, one can call the functions with reverse groups Y and X.
The treatment effect is then minus the estimated effect, and no adjustment is necessary for the standard error.

Generally,
`r waq(X,Y)$tau` and `r -waq(Y,X)$tau`
are not the same.
Under correct specification, they are asymptotically equivalent, but in finite samples they will differ based on which group is used to estimate the densities.
The same is also true when using `eif_additive`.

This package ports the Matlab code used in the simulations of Athey et al. (2021) to R.
Please use the Matlab replication code to exactly replicate the results and for implementations of other estimators studied by Athey et al. (2021).



## Reference

Susan Athey, Peter J. Bickel, Aiyou Chen, Guido W. Imbens, and Michael Pollmann. **Semiparametric Estimation of Treatment Effects in Randomized Experiments**. 2021. [[arxiv](https://arxiv.org/abs/2109.02603)]